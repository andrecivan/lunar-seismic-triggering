"""Efemerides lunares topocentricas vectorizadas (Skyfield + DE421).

NB09 hacia un loop de N llamadas individuales a Skyfield (~67000 evaluaciones
por iteracion MC -> minutos por iteracion). NB10 ya descubrio que Skyfield
acepta arrays. Aqui consolidamos ese patron como API publica.

Fisica:
- Posicion del observador: WGS84 (lat, lon, profundidad como elevacion negativa).
- Cuerpo objetivo: Luna (id JPL '301'), efemeride DE421.
- Salida: alt/az aparente y distancia geocentrica.

La marea cruda en CFS se modela como ~ 1/r^3 (gradiente del potencial).
"""
from __future__ import annotations
from typing import Iterable, Optional
from pathlib import Path
import numpy as np
from numpy.typing import NDArray
import pandas as pd

from skyfield.api import load, wgs84, Loader
from lunar_trigger.utils.paths import CACHE_DIR, NOTEBOOKS_DIR

# Cache compartido para el archivo de efemerides DE421 (~17 MB).
# Skyfield lo descarga una vez y luego lo lee de disco.
_EPHEM_DIR = CACHE_DIR / "skyfield"
_EPHEM_DIR.mkdir(parents=True, exist_ok=True)
_LOADER = Loader(str(_EPHEM_DIR))

_EPH = None
_TS = None


def _ensure_loaded() -> None:
    """Carga lazy de DE421 + escala de tiempos. Idempotente."""
    global _EPH, _TS
    if _EPH is None:
        # Re-usar de421.bsp si ya esta junto a los notebooks (compatibilidad con NB10).
        local_bsp = NOTEBOOKS_DIR / "de421.bsp"
        if local_bsp.exists() and not (_EPHEM_DIR / "de421.bsp").exists():
            import shutil
            shutil.copy(local_bsp, _EPHEM_DIR / "de421.bsp")
        _EPH = _LOADER("de421.bsp")
        _TS = _LOADER.timescale()


def moon_topocentric(
    latitudes: Iterable[float],
    longitudes: Iterable[float],
    depths_km: Iterable[float],
    times_utc: Iterable[pd.Timestamp],
) -> tuple[NDArray, NDArray, NDArray]:
    """Calcula alt/az/distancia de la Luna desde N observadores topocentricos.

    Vectorizada: una sola invocacion a Skyfield para todos los eventos.
    En benchmarks informales esto es ~50-100x mas rapido que llamar
    Skyfield N veces en un loop.

    Parameters
    ----------
    latitudes : array-like
        Latitudes geodeticas (grados).
    longitudes : array-like
        Longitudes (grados, este positivo).
    depths_km : array-like
        Profundidades del hipocentro (km, positivas hacia abajo).
        Se convierten internamente a elevacion negativa en metros.
    times_utc : array-like
        Timestamps UTC (datetime, pd.Timestamp o array de los mismos).

    Returns
    -------
    alt_deg, az_deg, dist_km : ndarray
        Altitud aparente, azimut y distancia topocentrica a la Luna.
    """
    _ensure_loaded()
    lats = np.asarray(latitudes, dtype=float)
    lons = np.asarray(longitudes, dtype=float)
    depths = np.asarray(depths_km, dtype=float)
    elevs_m = -depths * 1000.0  # depth km -> elevation -m

    earth = _EPH["earth"]
    moon = _EPH["moon"]

    # Skyfield acepta arrays de lat/lon -> crea un VectorObserver multi-punto.
    observer = earth + wgs84.latlon(lats, lons, elevation_m=elevs_m)

    # Normaliza tiempos a array de pd.Timestamp UTC (Skyfield exige tz-aware).
    dts = pd.to_datetime(list(times_utc))
    if dts.tz is None:
        dts = dts.tz_localize("UTC")
    else:
        dts = dts.tz_convert("UTC")
    t_arr = _TS.from_datetimes(dts.to_pydatetime())

    alt, az, dist = observer.at(t_arr).observe(moon).apparent().altaz()
    return alt.degrees, az.degrees, dist.km


def raw_tidal_components(
    alt_deg: NDArray,
    dist_km: NDArray,
) -> tuple[NDArray, NDArray]:
    """Componentes crudas (sigma_raw, tau_raw) de la marea ~ 1/r^3.

    Modelo simplificado: la amplitud decae como 1/r^3 (gradiente del potencial
    gravitacional), descompuesta en una componente vertical (sin alt) y una
    horizontal (cos alt). No incluye constantes fisicas porque solo importa
    la VARIACION temporal relativa para el test estadistico.
    """
    a_norm = 1.0 / (dist_km ** 3)
    alt_rad = np.radians(alt_deg)
    sigma_raw = a_norm * np.sin(alt_rad)
    tau_raw = a_norm * np.cos(alt_rad)
    return sigma_raw, tau_raw
