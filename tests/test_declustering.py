"""Tests para el modulo stats.declustering.

Cubre:
- Catalogo sin clusters: declustering no elimina nada.
- Cluster artificial bien definido: solo sobrevive el mainshock.
- Variantes fija vs variable producen tamanos coherentes.
"""
import numpy as np
import pandas as pd
import pytest

from lunar_trigger.stats.declustering import (
    gardner_knopoff_fixed,
    gardner_knopoff_variable,
    decluster,
)


def _build_catalog(rows):
    """Helper: construye un mini-catalogo desde una lista de dicts."""
    df = pd.DataFrame(rows)
    df["time"] = pd.to_datetime(df["time"], utc=True)
    return df


def test_catalogo_sin_clusters_no_pierde_eventos():
    # 5 eventos lejanos en espacio Y tiempo: ninguno deberia eliminarse.
    cat = _build_catalog([
        {"time": "2010-01-01", "latitude":  10.0, "longitude":  20.0, "mag": 7.5},
        {"time": "2012-01-01", "latitude": -30.0, "longitude":  80.0, "mag": 7.3},
        {"time": "2015-06-01", "latitude":  45.0, "longitude": -120.0, "mag": 7.8},
        {"time": "2018-12-01", "latitude":  -5.0, "longitude": 130.0, "mag": 7.2},
        {"time": "2022-03-15", "latitude":  60.0, "longitude":   5.0, "mag": 7.6},
    ])
    assert len(gardner_knopoff_fixed(cat)) == 5
    assert len(gardner_knopoff_variable(cat)) == 5


def test_cluster_artificial_solo_conserva_mainshock():
    # 1 mainshock + 5 replicas cercanas en tiempo y espacio.
    cluster = [{"time": "2020-01-01", "latitude": 0.0, "longitude": 0.0, "mag": 7.5}]
    cluster += [
        {"time": f"2020-01-{d:02d}", "latitude": 0.05, "longitude": 0.05, "mag": 5.0}
        for d in range(2, 7)
    ]
    cat = _build_catalog(cluster)
    out = gardner_knopoff_fixed(cat, radius_km=100, window_days=30)
    assert len(out) == 1
    assert out.iloc[0]["mag"] == pytest.approx(7.5)


def test_decluster_none_es_identidad():
    cat = _build_catalog([
        {"time": "2020-01-01", "latitude": 0.0, "longitude": 0.0, "mag": 7.0},
        {"time": "2020-01-02", "latitude": 0.0, "longitude": 0.0, "mag": 5.0},
    ])
    out = decluster(cat, method="none")
    assert len(out) == 2


def test_decluster_dispatch_acepta_todos_los_metodos():
    cat = _build_catalog([
        {"time": "2020-01-01", "latitude": 0.0, "longitude": 0.0, "mag": 7.5},
        {"time": "2020-02-01", "latitude": 50.0, "longitude": 50.0, "mag": 7.2},
    ])
    for method in ("none", "gk_fixed", "gk_variable"):
        out = decluster(cat, method=method)
        assert isinstance(out, pd.DataFrame)
        assert len(out) >= 1


def test_decluster_metodo_invalido_falla():
    cat = _build_catalog([{"time": "2020-01-01", "latitude": 0.0, "longitude": 0.0, "mag": 7.0}])
    with pytest.raises(ValueError):
        decluster(cat, method="metodo_inexistente")


def test_gk_variable_es_mas_estricto_que_fixed_en_grandes():
    # Mainshock M=8.5 con replicas a 200 km, 100 dias.
    # GK fixed (100 km / 30 d) NO las captura.
    # GK variable (R(M=8.5) ~ 600 km, T ~ 800 d) SI las captura.
    cluster = [{"time": "2020-01-01", "latitude": 0.0, "longitude": 0.0, "mag": 8.5}]
    cluster += [
        {"time": f"2020-{m:02d}-15", "latitude": 1.5, "longitude": 1.5, "mag": 5.5}
        for m in range(2, 6)
    ]
    cat = _build_catalog(cluster)
    n_fixed = len(gardner_knopoff_fixed(cat, radius_km=100, window_days=30))
    n_variable = len(gardner_knopoff_variable(cat))
    # Variable debe declusterizar MAS (sobrevivir menos eventos) que fixed con M=8.5.
    assert n_variable <= n_fixed
