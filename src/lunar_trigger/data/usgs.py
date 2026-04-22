"""Cliente USGS ComCat con sesion HTTP reutilizable y cache en disco.

Razones de existir:
- USGS limita peticiones por IP. Sin Session() abrimos un socket TCP por evento.
- Cualquier re-ejecucion del notebook bombardea la API innecesariamente.
- diskcache persiste resultados entre sesiones (sobrevive a reinicios del kernel).

API publica:
    fetch_catalog(...)  -> DataFrame del catalogo (parametrizado por fecha/magnitud).
    fetch_focal_mechanism(eventid)  -> dict con strike/dip/rake de NP1 y NP2 (o None).
"""
from __future__ import annotations
from typing import Optional
import time
import requests
import pandas as pd
import diskcache as dc

from lunar_trigger.utils.paths import CACHE_DIR

# Endpoints publicos de USGS ComCat (FDSN web services).
USGS_QUERY_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
USGS_EVENT_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={eid}&format=geojson"

# Cache compartido entre todas las funciones del modulo.
# Carpeta dedicada para no mezclarse con el cache de Skyfield.
_CACHE = dc.Cache(str(CACHE_DIR / "usgs"))

# Sesion HTTP unica (reusa conexion TCP entre llamadas, ~5x mas rapida).
_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "lunar-trigger/0.1 (research)"})


def _get_with_retries(url: str, params: Optional[dict] = None,
                      max_retries: int = 3, backoff: float = 1.5) -> requests.Response:
    """GET resiliente con backoff exponencial.

    USGS responde 503 ocasionalmente bajo carga; reintentar es casi siempre suficiente.
    """
    last_exc = None
    for attempt in range(max_retries):
        try:
            r = _SESSION.get(url, params=params, timeout=30)
            r.raise_for_status()
            return r
        except (requests.RequestException, requests.HTTPError) as e:
            last_exc = e
            if attempt < max_retries - 1:
                time.sleep(backoff ** attempt)
    raise RuntimeError(f"USGS GET fallo tras {max_retries} intentos: {last_exc}")


def fetch_catalog(
    starttime: str,
    endtime: str,
    minmagnitude: float = 7.0,
    maxdepth: Optional[float] = 70.0,
    orderby: str = "time",
    limit: int = 20000,
) -> pd.DataFrame:
    """Descarga un catalogo USGS y lo devuelve como DataFrame con columnas estandar.

    El resultado se cachea en disco usando los parametros como clave: si vuelves
    a llamar con los mismos argumentos, no se hace request a USGS.

    Parameters
    ----------
    starttime, endtime : str
        Fechas ISO (ej: "1995-01-01", "2024-12-31").
    minmagnitude : float
        Magnitud minima (default 7.0 para mega-terremotos).
    maxdepth : float | None
        Profundidad maxima en km. Default 70 (sismos crustales/superficiales).
    """
    cache_key = ("catalog", starttime, endtime, minmagnitude, maxdepth, orderby, limit)
    if cache_key in _CACHE:
        return _CACHE[cache_key].copy()

    params = {
        "format": "csv",
        "starttime": starttime,
        "endtime": endtime,
        "minmagnitude": minmagnitude,
        "orderby": orderby,
        "limit": limit,
    }
    if maxdepth is not None:
        params["maxdepth"] = maxdepth

    r = _get_with_retries(USGS_QUERY_URL, params=params)

    from io import StringIO
    df = pd.read_csv(StringIO(r.text))
    # Normaliza el campo de tiempo (USGS mezcla microsegundos en versiones recientes).
    if "time" in df.columns:
        df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")

    _CACHE[cache_key] = df.copy()
    return df


def fetch_focal_mechanism(event_id: str) -> Optional[dict]:
    """Recupera el tensor momento (strike/dip/rake de NP1 y NP2) para un evento.

    USGS guarda los mecanismos focales como sub-productos del evento principal.
    Devuelve None si el evento no tiene mecanismo publicado (comun en sismos
    pequenos o muy recientes).

    Returns
    -------
    dict | None
        {strike1, dip1, rake1, strike2, dip2, rake2} o None.
    """
    cache_key = ("focal", event_id)
    if cache_key in _CACHE:
        return _CACHE[cache_key]

    url = USGS_EVENT_URL.format(eid=event_id)
    try:
        r = _get_with_retries(url)
    except RuntimeError:
        _CACHE[cache_key] = None
        return None

    data = r.json()
    products = data.get("properties", {}).get("products", {})
    mts = products.get("moment-tensor") or products.get("focal-mechanism")
    if not mts:
        _CACHE[cache_key] = None
        return None

    # Tomamos el primer producto (preferred). USGS ya lo ordena por preferencia.
    props = mts[0].get("properties", {})
    try:
        result = {
            "strike1": float(props["nodal-plane-1-strike"]),
            "dip1":    float(props["nodal-plane-1-dip"]),
            "rake1":   float(props["nodal-plane-1-rake"]),
            "strike2": float(props["nodal-plane-2-strike"]),
            "dip2":    float(props["nodal-plane-2-dip"]),
            "rake2":   float(props["nodal-plane-2-rake"]),
        }
    except (KeyError, ValueError, TypeError):
        result = None

    _CACHE[cache_key] = result
    return result


def clear_cache() -> None:
    """Vacia el cache USGS. Util si las APIs cambiaron o necesitas re-descargar."""
    _CACHE.clear()
