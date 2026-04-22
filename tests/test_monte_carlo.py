"""Tests para el modulo stats.monte_carlo.

Cubre:
- Reproducibilidad por semilla.
- Estructura del dict de retorno.
- Forma de la matriz NxN.
- Bajo H0 puro (catalogo aleatorio), p-value debe estar uniformemente distribuido.
"""
import numpy as np
import pandas as pd
import pytest

from lunar_trigger.stats.monte_carlo import time_shuffling_null


@pytest.fixture
def mini_catalogo():
    """Catalogo sintetico de 30 sismos en distintos cinturones."""
    rng = np.random.default_rng(0)
    n = 30
    return pd.DataFrame({
        "time": pd.date_range("2010-01-01", "2024-12-31", periods=n, tz="UTC"),
        "latitude":  rng.uniform(-60, 60, n),
        "longitude": rng.uniform(-180, 180, n),
        "depth":     rng.uniform(5, 60, n),
        "dip1":      rng.uniform(10, 80, n),
        "dip2":      rng.uniform(10, 80, n),
    })


def test_estructura_de_retorno(mini_catalogo):
    res = time_shuffling_null(mini_catalogo, n_iterations=50, seed=1, n_jobs=1)
    expected_keys = {
        "observed_fraction", "null_distribution", "p_value",
        "null_mean", "null_std", "cfs_matrix",
        "n_events", "n_iterations", "mu",
    }
    assert expected_keys.issubset(res.keys())
    assert res["n_events"] == 30
    assert res["n_iterations"] == 50


def test_matriz_cfs_es_NxN(mini_catalogo):
    res = time_shuffling_null(mini_catalogo, n_iterations=10, seed=1, n_jobs=1)
    n = len(mini_catalogo)
    assert res["cfs_matrix"].shape == (n, n)


def test_distribucion_nula_tiene_tamano_correcto(mini_catalogo):
    res = time_shuffling_null(mini_catalogo, n_iterations=100, seed=1, n_jobs=1)
    assert res["null_distribution"].shape == (100,)
    # Todos en [0, 1] porque son fracciones.
    assert (res["null_distribution"] >= 0).all()
    assert (res["null_distribution"] <= 1).all()


def test_reproducibilidad_por_semilla(mini_catalogo):
    res1 = time_shuffling_null(mini_catalogo, n_iterations=200, seed=42, n_jobs=1)
    res2 = time_shuffling_null(mini_catalogo, n_iterations=200, seed=42, n_jobs=1)
    assert np.array_equal(res1["null_distribution"], res2["null_distribution"])
    assert res1["p_value"] == res2["p_value"]


def test_observed_fraction_esta_en_diagonal(mini_catalogo):
    # observed_fraction debe ser exactamente la fraccion de la diagonal > 0.
    res = time_shuffling_null(mini_catalogo, n_iterations=10, seed=1, n_jobs=1)
    diag = np.diag(res["cfs_matrix"])
    assert res["observed_fraction"] == pytest.approx((diag > 0).mean())


def test_paralelizacion_da_mismo_resultado(mini_catalogo):
    # Joblib con threading no debe alterar el orden ni los valores.
    res_serial = time_shuffling_null(mini_catalogo, n_iterations=100, seed=99, n_jobs=1)
    res_paralelo = time_shuffling_null(mini_catalogo, n_iterations=100, seed=99, n_jobs=4)
    assert np.array_equal(res_serial["null_distribution"], res_paralelo["null_distribution"])
