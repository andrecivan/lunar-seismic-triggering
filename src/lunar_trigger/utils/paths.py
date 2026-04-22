"""Rutas absolutas del proyecto.

Resuelve la raiz del repositorio a partir de la ubicacion de este archivo,
de modo que cualquier notebook (cualquiera que sea su CWD) pueda importar
estas constantes y leer/escribir datos sin romperse por cambios de directorio.
"""
from __future__ import annotations
from pathlib import Path

# src/lunar_trigger/utils/paths.py -> raiz = 4 niveles arriba
PROJECT_ROOT: Path = Path(__file__).resolve().parents[3]

# Datos crudos: descargados directamente de fuentes externas (USGS, JPL).
DATA_RAW: Path = PROJECT_ROOT / "data" / "raw"

# Datos derivados: salidas intermedias de los notebooks (con efemerides, CFS, etc).
DATA_PROCESSED: Path = PROJECT_ROOT / "data" / "processed"

# Cache en disco (USGS, Skyfield) para evitar llamadas repetidas a redes externas.
CACHE_DIR: Path = PROJECT_ROOT / "data" / "cache"

# Resultados de publicacion: figuras, mapas, tablas estadisticas.
RESULTS_FIGURES: Path = PROJECT_ROOT / "results" / "figures"
RESULTS_MAPS: Path = PROJECT_ROOT / "results" / "maps"
RESULTS_STATS: Path = PROJECT_ROOT / "results" / "statistics"

# Notebooks (a veces se necesita ubicar el .bsp de Skyfield aqui).
NOTEBOOKS_DIR: Path = PROJECT_ROOT / "notebooks"


def ensure_dirs() -> None:
    """Crea todos los directorios estandar si no existen.

    Idempotente: seguro de llamar multiples veces. Util para inicializar
    un entorno fresco antes de ejecutar el pipeline completo.
    """
    for d in (
        DATA_RAW,
        DATA_PROCESSED,
        CACHE_DIR,
        RESULTS_FIGURES,
        RESULTS_MAPS,
        RESULTS_STATS,
    ):
        d.mkdir(parents=True, exist_ok=True)
