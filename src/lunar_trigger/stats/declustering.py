"""Declustering de catalogos sismicos.

Por que importa: el modelo MC asume independencia temporal de los eventos.
Si dos sismos cercanos en tiempo y espacio son mainshock + replica, comparten
una fisica COMUN (transferencia de Coulomb estatico) que NO depende de la marea.
Sin declustering, el catalogo viola la hipotesis nula y produce p-values
inflados.

Implementaciones:
- gardner_knopoff_fixed:    ventanas (R, T) constantes (rapido, conservador).
- gardner_knopoff_variable: R(M), T(M) crecen con la magnitud (Gardner-Knopoff
  1974 + Helmstetter & Sornette parametrization). Estandar de la literatura.

NB10 implementaba ventanas fijas con dos loops anidados O(N^2). Aqui usamos
scipy.spatial.cKDTree para vectorizar la busqueda de vecinos espaciales.
"""
from __future__ import annotations
from typing import Literal
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# Radio terrestre medio (km), suficiente para la aproximacion plana local
# en distancias < 500 km (error < 1%).
EARTH_RADIUS_KM = 6371.0


def _haversine_km(lat1, lon1, lat2, lon2):
    """Distancia great-circle entre dos puntos en km."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def _to_xyz(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    """Convierte lat/lon a coordenadas cartesianas en una esfera de R=1 (km after scale).

    cKDTree mide distancia euclidea; sobre la esfera unitaria, la cuerda 3D
    aproxima muy bien el arco para distancias pequenas (< 1000 km, error < 1%).
    Multiplicamos por R para devolver km de cuerda.
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    x = np.cos(lat) * np.cos(lon) * EARTH_RADIUS_KM
    y = np.cos(lat) * np.sin(lon) * EARTH_RADIUS_KM
    z = np.sin(lat) * EARTH_RADIUS_KM
    return np.column_stack([x, y, z])


def gardner_knopoff_fixed(
    catalog: pd.DataFrame,
    radius_km: float = 100.0,
    window_days: float = 30.0,
    time_col: str = "time",
    lat_col: str = "latitude",
    lon_col: str = "longitude",
    mag_col: str = "mag",
) -> pd.DataFrame:
    """Declustering con ventanas (R, T) fijas. Replica de NB10 pero vectorizado.

    Algoritmo:
    1. Ordena por magnitud descendente (mainshocks primero).
    2. Para cada mainshock candidato, marca como replica todo sismo posterior
       a el dentro de (radius_km, window_days).
    3. Devuelve solo los mainshocks no marcados.

    Esto es O(N log N) gracias a cKDTree (vs O(N^2) del loop original).

    Returns
    -------
    DataFrame
        Subconjunto del catalogo con solo mainshocks (declusterizado).
    """
    df = catalog.copy().reset_index(drop=True)
    df[time_col] = pd.to_datetime(df[time_col], utc=True)

    # Orden de procesamiento: magnitud descendente (mainshocks dominan).
    order = df.sort_values(mag_col, ascending=False).index.to_numpy()

    xyz = _to_xyz(df[lat_col].to_numpy(), df[lon_col].to_numpy())
    tree = cKDTree(xyz)
    times = df[time_col].to_numpy()

    is_mainshock = np.ones(len(df), dtype=bool)
    window_td = pd.Timedelta(days=window_days).to_timedelta64()

    for i in order:
        if not is_mainshock[i]:
            continue
        # Vecinos espaciales dentro de radius_km (cuerda 3D ~ arco para r<1000km).
        neighbors = tree.query_ball_point(xyz[i], r=radius_km)
        for j in neighbors:
            if j == i or not is_mainshock[j]:
                continue
            # Solo eventos POSTERIORES y dentro de la ventana temporal.
            dt = times[j] - times[i]
            if dt > np.timedelta64(0) and dt <= window_td:
                is_mainshock[j] = False

    return df.loc[is_mainshock].reset_index(drop=True)


def _gk_window_helmstetter(mag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Ventanas R(M), T(M) parametrizadas por Helmstetter & Sornette (2003)
    sobre Gardner-Knopoff (1974).

    R en km, T en dias.
    """
    radius = 10 ** (0.1238 * mag + 0.983)
    days = np.where(
        mag >= 6.5,
        10 ** (0.032 * mag + 2.7389),
        10 ** (0.5409 * mag - 0.547),
    )
    return radius, days


def gardner_knopoff_variable(
    catalog: pd.DataFrame,
    time_col: str = "time",
    lat_col: str = "latitude",
    lon_col: str = "longitude",
    mag_col: str = "mag",
) -> pd.DataFrame:
    """Declustering con ventanas R(M), T(M) variables segun magnitud.

    Estandar de la literatura sismologica (Gardner-Knopoff 1974, parametros
    actualizados por Helmstetter & Sornette 2003). Mas estricto que ventanas
    fijas para mainshocks grandes (un M=8.5 declusteriza ~1000 km / ~3 anos).

    El ancho de ventana se toma del MAINSHOCK candidato, no del sismo a evaluar.
    """
    df = catalog.copy().reset_index(drop=True)
    df[time_col] = pd.to_datetime(df[time_col], utc=True)

    mags = df[mag_col].to_numpy()
    order = df.sort_values(mag_col, ascending=False).index.to_numpy()

    radii_km, days = _gk_window_helmstetter(mags)
    windows_td = pd.to_timedelta(days, unit="D").to_numpy()

    xyz = _to_xyz(df[lat_col].to_numpy(), df[lon_col].to_numpy())
    tree = cKDTree(xyz)
    times = df[time_col].to_numpy()
    is_mainshock = np.ones(len(df), dtype=bool)

    for i in order:
        if not is_mainshock[i]:
            continue
        r_i = radii_km[i]
        w_i = windows_td[i]
        neighbors = tree.query_ball_point(xyz[i], r=r_i)
        for j in neighbors:
            if j == i or not is_mainshock[j]:
                continue
            dt = times[j] - times[i]
            # Replica si esta DESPUES del mainshock dentro de la ventana,
            # O foreshock dentro de una ventana mas corta (1/4 del total).
            if dt > np.timedelta64(0) and dt <= w_i:
                is_mainshock[j] = False
            elif dt < np.timedelta64(0) and -dt <= (w_i / 4):
                is_mainshock[j] = False

    return df.loc[is_mainshock].reset_index(drop=True)


def decluster(
    catalog: pd.DataFrame,
    method: Literal["none", "gk_fixed", "gk_variable"] = "gk_variable",
    **kwargs,
) -> pd.DataFrame:
    """Despachador comodo. method='none' devuelve el catalogo intacto."""
    if method == "none":
        return catalog.copy()
    if method == "gk_fixed":
        return gardner_knopoff_fixed(catalog, **kwargs)
    if method == "gk_variable":
        return gardner_knopoff_variable(catalog, **kwargs)
    raise ValueError(f"Metodo desconocido: {method}")
