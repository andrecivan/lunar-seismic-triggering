"""Monte Carlo time-shuffling con pre-computacion matricial.

Idea clave para acelerar:
- Lo unico que cambia entre iteraciones MC es la asignacion tiempo<->ubicacion.
- Las ubicaciones (lat, lon, depth) y los planos nodales (dip1, dip2) son fijos.
- Si pre-computamos una matriz NxN de "que pasaria si la ubicacion i hubiera
  ocurrido en el tiempo j", cada iteracion MC se reduce a una indexacion por
  permutacion + agregacion. Sin re-llamar a Skyfield.

NB10 ya vectorizo Skyfield POR iteracion (1 llamada con array de tiempos).
Esto va un paso mas alla: 1 sola llamada con N*N pares (i_loc, j_time),
y luego cada iteracion es solo numpy fancy indexing.

Tradeoff:
- Memoria: NxN floats. Para N=335 son ~900 KB. Trivial.
- Tiempo Skyfield: 1 llamada con N^2 evaluaciones. Para N=335 son 112225
  evaluaciones, similar al costo de NB10 con 1000 iteraciones x 335 (=335000),
  pero amortizado en una sola llamada (~3x mas rapido por TCP+overhead).
- Total: ~10x mas rapido que NB10 cuando se combina con pre-computacion +
  joblib en paralelo.
"""
from __future__ import annotations
from typing import Callable, Optional
import numpy as np
from numpy.typing import NDArray
import pandas as pd
from joblib import Parallel, delayed

from lunar_trigger.data.ephemeris import moon_topocentric, raw_tidal_components
from lunar_trigger.physics.coulomb import calc_cfs_unbiased, DEFAULT_FRICTION


def _build_cfs_matrix(
    latitudes: NDArray,
    longitudes: NDArray,
    depths_km: NDArray,
    dip1_deg: NDArray,
    dip2_deg: NDArray,
    times_utc: pd.DatetimeIndex,
    mu: float,
) -> NDArray:
    """Pre-computa la matriz NxN de Delta CFS para todos los pares (loc_i, time_j).

    Returns
    -------
    M : ndarray shape (N, N)
        M[i, j] = Delta CFS si el evento en ubicacion i hubiera ocurrido en tiempo j.
    """
    n = len(latitudes)

    # Producto cartesiano (loc_i, time_j) aplanado en (N*N,) para una sola llamada Skyfield.
    lat_grid = np.repeat(latitudes, n)            # [l0, l0, ..., l1, l1, ..., l_{n-1}, ...]
    lon_grid = np.repeat(longitudes, n)
    dep_grid = np.repeat(depths_km, n)
    dip1_grid = np.repeat(dip1_deg, n)
    dip2_grid = np.repeat(dip2_deg, n)
    time_grid = np.tile(times_utc, n)             # [t0, t1, ..., t_{n-1}, t0, t1, ...]

    alt, _, dist = moon_topocentric(lat_grid, lon_grid, dep_grid, time_grid)
    sigma, tau = raw_tidal_components(alt, dist)
    cfs_flat = calc_cfs_unbiased(sigma, tau, dip1_grid, dip2_grid, mu=mu)
    return cfs_flat.reshape(n, n)


def _iter_fraction_favorable(
    cfs_matrix: NDArray,
    perm: NDArray,
) -> float:
    """Una iteracion MC: aplica permutacion de tiempos y mide % CFS > 0.

    perm[i] = indice de tiempo asignado a la ubicacion i.
    Diagonal de M[arange(N), perm] = CFS reasignado para esta iteracion.
    """
    n = cfs_matrix.shape[0]
    cfs_iter = cfs_matrix[np.arange(n), perm]
    return float(np.mean(cfs_iter > 0.0))


def time_shuffling_null(
    catalog: pd.DataFrame,
    n_iterations: int = 1000,
    mu: float = DEFAULT_FRICTION,
    seed: Optional[int] = 42,
    n_jobs: int = -1,
    columns: Optional[dict] = None,
) -> dict:
    """Modelo nulo Monte Carlo barajando los tiempos de ocurrencia.

    Pipeline:
    1. Pre-computa la matriz NxN de CFS (una unica llamada vectorizada a Skyfield).
    2. Para cada iteracion: genera una permutacion aleatoria de los N tiempos
       y mide la fraccion de eventos con CFS > 0 (favorables).
    3. Devuelve la distribucion nula completa (n_iterations valores).

    Parameters
    ----------
    catalog : DataFrame
        Catalogo declusterizado con columnas: time, latitude, longitude, depth,
        dip1, dip2 (configurables via `columns`).
    n_iterations : int
        Numero de permutaciones MC (1000 -> resolucion p~0.001).
    mu : float
        Coeficiente de friccion para el calculo CFS.
    seed : int | None
        Semilla del PRNG para reproducibilidad (None = aleatorio puro).
    n_jobs : int
        Procesos paralelos (-1 = todos los cores). joblib threading suele bastar
        porque la operacion por iteracion es 100% numpy (libera el GIL).
    columns : dict | None
        Mapeo de nombres si el catalogo usa nombres distintos (ej: {'depth': 'depth_km'}).

    Returns
    -------
    dict con:
        observed_fraction : fraccion CFS>0 de la asignacion REAL (sin permutar).
        null_distribution : np.ndarray de fracciones bajo H0.
        p_value           : fraccion del nulo >= observado (one-sided).
        null_mean, null_std : estadisticos del nulo.
        cfs_matrix        : matriz NxN (util para diagnosticos).
    """
    cols = {"time": "time", "latitude": "latitude", "longitude": "longitude",
            "depth": "depth", "dip1": "dip1", "dip2": "dip2"}
    if columns:
        cols.update(columns)

    df = catalog.copy()
    times = pd.to_datetime(df[cols["time"]], utc=True)
    lats = df[cols["latitude"]].to_numpy(dtype=float)
    lons = df[cols["longitude"]].to_numpy(dtype=float)
    deps = df[cols["depth"]].to_numpy(dtype=float)
    dip1 = df[cols["dip1"]].to_numpy(dtype=float)
    dip2 = df[cols["dip2"]].to_numpy(dtype=float)
    n = len(df)

    # Paso 1: matriz NxN de CFS (una sola llamada Skyfield).
    M = _build_cfs_matrix(lats, lons, deps, dip1, dip2,
                          pd.DatetimeIndex(times), mu)

    # Observado real = diagonal de M (cada loc i con su tiempo i original).
    observed = float(np.mean(np.diag(M) > 0.0))

    # Paso 2: distribucion nula via permutaciones.
    rng = np.random.default_rng(seed)
    # Pre-genera todas las permutaciones (compacto, evita race conditions con joblib).
    perms = [rng.permutation(n) for _ in range(n_iterations)]

    # joblib con threading: cada iteracion es numpy puro, libera GIL.
    null_dist = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_iter_fraction_favorable)(M, p) for p in perms
    )
    null_dist = np.array(null_dist, dtype=float)

    # P-value empirico one-sided: fraccion del nulo >= observado.
    p_value = float(np.mean(null_dist >= observed))

    return {
        "observed_fraction": observed,
        "null_distribution": null_dist,
        "p_value": p_value,
        "null_mean": float(null_dist.mean()),
        "null_std": float(null_dist.std(ddof=1)),
        "cfs_matrix": M,
        "n_events": n,
        "n_iterations": n_iterations,
        "mu": mu,
    }
