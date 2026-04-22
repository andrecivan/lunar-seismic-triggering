"""Lunar Tidal Triggering of Earthquakes.

Paquete que centraliza el pipeline de la investigacion:
- data: ingesta de catalogos sismicos (USGS) y efemerides lunares (Skyfield/JPL).
- physics: modelos geomecanicos (Criterio de Falla de Coulomb).
- stats: validacion estadistica (Monte Carlo, declustering, tests circulares).
- utils: utilidades transversales (paths, logging, IO).

Filosofia: cada notebook se vuelve un script delgado que importa funciones puras
testeadas. La fisica vive en el paquete, la narrativa vive en los notebooks.
"""

__version__ = "0.1.0"
__author__ = "Ivan Andres Mena Contreras"
