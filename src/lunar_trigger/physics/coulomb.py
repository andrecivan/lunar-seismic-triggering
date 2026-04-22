"""Criterio de Falla de Coulomb proyectado a planos nodales reales.

Centraliza el calculo de Delta CFS = Delta tau + mu * Delta sigma_n
para resolver la ambiguedad nodal del mecanismo focal de forma desesgada
(promedio de los dos planos ortogonales).

Esta funcion es el corazon fisico de la investigacion: NB09 y NB10 dependen
de ella. Antes existia duplicada en ambos notebooks; cualquier divergencia
producia resultados estadisticos inconsistentes. Ahora vive en un solo lugar.
"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike, NDArray

# Coeficiente de friccion estatica por defecto.
# 0.4 es valor estandar de la literatura para fallas de subduccion (Cocco & Rice 2002).
DEFAULT_FRICTION: float = 0.4


def project_to_fault(
    sigma_raw: ArrayLike,
    tau_raw: ArrayLike,
    dip_deg: ArrayLike,
) -> tuple[NDArray, NDArray]:
    """Proyecta el vector de marea topocentrico sobre un plano de falla.

    Aplica una rotacion trigonometrica 2D usando el angulo de buzamiento (dip)
    para descomponer el esfuerzo crudo en sus componentes normal (sigma_n) y
    cortante (tau) RESPECTO al plano de la falla.

    Parameters
    ----------
    sigma_raw : array_like
        Componente vertical de la perturbacion gravitacional lunar (proporcional a sin(alt)/r^3).
    tau_raw : array_like
        Componente horizontal (proporcional a cos(alt)/r^3).
    dip_deg : array_like
        Angulo de buzamiento del plano de falla, en grados.

    Returns
    -------
    fault_normal, fault_shear : ndarray
        Esfuerzo normal (positivo = "unclamping") y cortante sobre la falla.
    """
    dip_rad = np.radians(np.asarray(dip_deg, dtype=float))
    s = np.asarray(sigma_raw, dtype=float)
    t = np.asarray(tau_raw, dtype=float)
    fault_normal = s * np.cos(dip_rad) - t * np.sin(dip_rad)
    fault_shear = s * np.sin(dip_rad) + t * np.cos(dip_rad)
    return fault_normal, fault_shear


def cfs_single_plane(
    sigma_raw: ArrayLike,
    tau_raw: ArrayLike,
    dip_deg: ArrayLike,
    mu: float = DEFAULT_FRICTION,
) -> NDArray:
    """Delta CFS proyectado sobre UN unico plano nodal.

    Coulomb Failure Stress = tau + mu * sigma_n
    (signo positivo = el plano se acerca a su umbral de ruptura).
    """
    fault_normal, fault_shear = project_to_fault(sigma_raw, tau_raw, dip_deg)
    return fault_shear + mu * fault_normal


def calc_cfs_unbiased(
    sigma_raw: ArrayLike,
    tau_raw: ArrayLike,
    dip1_deg: ArrayLike,
    dip2_deg: ArrayLike,
    mu: float = DEFAULT_FRICTION,
) -> NDArray:
    """Delta CFS DESESGADO promediando los dos planos nodales ortogonales.

    Cada mecanismo focal admite dos planos matematicamente equivalentes (NP1, NP2).
    Sin informacion adicional (afterslip, geologia local) NO podemos elegir.

    - Tomar max(CFS1, CFS2) introduce SESGO POSITIVO: como la marea oscila,
      siempre habra un plano "favorable" por azar. Inflaria el % de eventos
      favorecidos artificialmente.
    - Promediar (CFS1 + CFS2) / 2 asume ignorancia total y es la opcion conservadora.

    Esta es la senal "real" del 64.86% reportada en el manuscrito GRL (Fase 4).

    Parameters
    ----------
    sigma_raw, tau_raw : array_like
        Componentes crudas de la marea topocentrica (1/r^3).
    dip1_deg, dip2_deg : array_like
        Buzamientos de los dos planos nodales (NP1 y NP2) en grados.
    mu : float
        Coeficiente de friccion estatica (default 0.4).

    Returns
    -------
    ndarray
        Delta CFS promediado, una entrada por evento.
    """
    cfs1 = cfs_single_plane(sigma_raw, tau_raw, dip1_deg, mu=mu)
    cfs2 = cfs_single_plane(sigma_raw, tau_raw, dip2_deg, mu=mu)
    return (cfs1 + cfs2) / 2.0


def calc_cfs_biased_max(
    sigma_raw: ArrayLike,
    tau_raw: ArrayLike,
    dip1_deg: ArrayLike,
    dip2_deg: ArrayLike,
    mu: float = DEFAULT_FRICTION,
) -> NDArray:
    """Variante SESGADA: max(CFS1, CFS2). Usada solo en NB08 con fines historicos.

    Mantiene paridad con la formulacion original del notebook 08, pero NO debe
    usarse para conclusiones: introduce sesgo de seleccion positivo. Documentada
    aqui para que NB08 pueda importarla y reproducir su resultado preliminar.
    """
    cfs1 = cfs_single_plane(sigma_raw, tau_raw, dip1_deg, mu=mu)
    cfs2 = cfs_single_plane(sigma_raw, tau_raw, dip2_deg, mu=mu)
    return np.maximum(cfs1, cfs2)
