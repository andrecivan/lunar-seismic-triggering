"""Tests para el modulo physics.coulomb.

Cubre:
- Comportamientos triviales (dip=0, dip=90, mu=0).
- Equivalencia entre formulacion vectorial nueva y la formula original duplicada.
- Conmutatividad del promedio nodal: NP1+NP2 == NP2+NP1.
"""
import numpy as np
import pytest

from lunar_trigger.physics.coulomb import (
    project_to_fault,
    cfs_single_plane,
    calc_cfs_unbiased,
    calc_cfs_biased_max,
    DEFAULT_FRICTION,
)


def _legacy_cfs_unbiased(sigma_raw, tau_raw, dip1, dip2, mu=0.4):
    """Replica EXACTA de la funcion duplicada en NB09 y NB10. Usada como oraculo."""
    dip1_rad = np.radians(dip1); dip2_rad = np.radians(dip2)
    n1 = sigma_raw * np.cos(dip1_rad) - tau_raw * np.sin(dip1_rad)
    s1 = sigma_raw * np.sin(dip1_rad) + tau_raw * np.cos(dip1_rad)
    cfs1 = s1 + (mu * n1)
    n2 = sigma_raw * np.cos(dip2_rad) - tau_raw * np.sin(dip2_rad)
    s2 = sigma_raw * np.sin(dip2_rad) + tau_raw * np.cos(dip2_rad)
    cfs2 = s2 + (mu * n2)
    return (cfs1 + cfs2) / 2.0


def test_project_dip_zero_es_identidad():
    # Dip=0 -> el plano de falla es horizontal: sigma queda como esta, tau queda como esta.
    sigma = np.array([1.0, 2.0, -3.0])
    tau = np.array([0.5, -1.0, 4.0])
    n, s = project_to_fault(sigma, tau, dip_deg=0.0)
    assert np.allclose(n, sigma)
    assert np.allclose(s, tau)


def test_project_dip_noventa_intercambia():
    # Dip=90 -> rotacion de 90 grados: las componentes se intercambian con signo.
    sigma = np.array([1.0])
    tau = np.array([2.0])
    n, s = project_to_fault(sigma, tau, dip_deg=90.0)
    # cos(90)=0, sin(90)=1
    # n = sigma*0 - tau*1 = -tau
    # s = sigma*1 + tau*0 = sigma
    assert np.allclose(n, -tau)
    assert np.allclose(s, sigma)


def test_cfs_single_friccion_cero_devuelve_solo_shear():
    # mu=0 -> CFS = tau pura (sin contribucion normal).
    sigma = np.array([5.0])
    tau = np.array([3.0])
    _, shear_real = project_to_fault(sigma, tau, dip_deg=30.0)
    cfs = cfs_single_plane(sigma, tau, dip_deg=30.0, mu=0.0)
    assert np.allclose(cfs, shear_real)


def test_unbiased_paridad_con_legacy_aleatorio():
    # Genera 10000 casos aleatorios, compara nueva implementacion vs duplicada original.
    rng = np.random.default_rng(123)
    n = 10000
    sigma = rng.normal(0, 1e-6, n)
    tau = rng.normal(0, 1e-6, n)
    d1 = rng.uniform(5, 85, n)
    d2 = rng.uniform(5, 85, n)
    nueva = calc_cfs_unbiased(sigma, tau, d1, d2, mu=DEFAULT_FRICTION)
    legacy = _legacy_cfs_unbiased(sigma, tau, d1, d2, mu=DEFAULT_FRICTION)
    assert np.allclose(nueva, legacy, atol=1e-20)


def test_unbiased_es_conmutativo_en_planos():
    # NP1+NP2 = NP2+NP1 (promedio simetrico).
    rng = np.random.default_rng(7)
    n = 500
    sigma = rng.normal(0, 1, n)
    tau = rng.normal(0, 1, n)
    d1 = rng.uniform(5, 85, n)
    d2 = rng.uniform(5, 85, n)
    a = calc_cfs_unbiased(sigma, tau, d1, d2)
    b = calc_cfs_unbiased(sigma, tau, d2, d1)
    assert np.allclose(a, b)


def test_biased_max_es_siempre_geq_unbiased():
    # max(CFS1, CFS2) >= mean(CFS1, CFS2) por definicion.
    rng = np.random.default_rng(0)
    sigma = rng.normal(0, 1, 1000)
    tau = rng.normal(0, 1, 1000)
    d1 = rng.uniform(5, 85, 1000)
    d2 = rng.uniform(5, 85, 1000)
    biased = calc_cfs_biased_max(sigma, tau, d1, d2)
    unbiased = calc_cfs_unbiased(sigma, tau, d1, d2)
    assert (biased >= unbiased - 1e-12).all()


@pytest.mark.parametrize("mu", [0.0, 0.2, 0.4, 0.6, 0.8])
def test_friccion_modula_linealmente(mu):
    # CFS = tau + mu*sigma_n, debe ser lineal en mu para cualquier evento.
    sigma = np.array([1.5])
    tau = np.array([0.5])
    n_proy, s_proy = project_to_fault(sigma, tau, dip_deg=45.0)
    cfs = cfs_single_plane(sigma, tau, dip_deg=45.0, mu=mu)
    assert np.allclose(cfs, s_proy + mu * n_proy)
