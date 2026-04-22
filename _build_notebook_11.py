"""Constructor del Notebook 11: matriz de sensibilidad de declustering.

Genera notebooks/11_declustering_sensitivity.ipynb con celdas que ejercitan
el paquete `lunar_trigger` recien creado: comparan p-value bajo distintos
metodos de declustering para evaluar robustez del resultado.

Script de uso unico (regenerar el .ipynb si cambia el contenido).
"""
from __future__ import annotations
import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()

cells = []

cells.append(nbf.v4.new_markdown_cell("""\
# Notebook 11 - Matriz de Sensibilidad: Declustering vs Resultado MC

**Objetivo cientifico:** evaluar la robustez del resultado del manuscrito GRL
(Notebook 10) frente a la eleccion del metodo de declustering. Si la senal
del *Lunar Tidal Triggering* es real, el `p-value` debe sobrevivir a multiples
formas razonables de eliminar replicas.

**Metodologia:**
1. Catalogo base: `data/processed/earthquakes_global_robust.csv` (333 sismos M>=7).
2. Tres niveles de declustering: ninguno, Gardner-Knopoff fijo (100 km / 30 d),
   Gardner-Knopoff variable (Helmstetter & Sornette 2003).
3. Para cada catalogo declusterizado, Monte Carlo time-shuffling 1000x con la
   funcion centralizada `lunar_trigger.stats.monte_carlo.time_shuffling_null`.
4. Tabla resumen de N, fraccion observada, media nula, p-value.

**Esto reemplaza la duplicacion de codigo de NB09/NB10**: ahora todo importa
desde el paquete instalable `lunar_trigger` (instalado con `pip install -e .`).
"""))

cells.append(nbf.v4.new_code_cell("""\
# Imports del paquete y dependencias estandar
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lunar_trigger.utils.paths import DATA_PROCESSED, RESULTS_FIGURES, RESULTS_STATS, ensure_dirs
from lunar_trigger.stats.declustering import gardner_knopoff_fixed, gardner_knopoff_variable
from lunar_trigger.stats.monte_carlo import time_shuffling_null
from lunar_trigger.physics.coulomb import DEFAULT_FRICTION

ensure_dirs()
print('Friccion estatica por defecto: mu =', DEFAULT_FRICTION)\
"""))

cells.append(nbf.v4.new_markdown_cell("""\
## 1. Carga del catalogo base

`earthquakes_global_robust.csv` ya contiene los dos planos nodales (`dip1`, `dip2`)
para el calculo desesgado de Coulomb stress.
"""))

cells.append(nbf.v4.new_code_cell("""\
# Carga catalogo con planos nodales completos
df_base = pd.read_csv(DATA_PROCESSED / 'earthquakes_global_robust.csv')
df_base['time'] = pd.to_datetime(df_base['time'], utc=True, format='ISO8601')

# Renombra 'magnitude' -> 'mag' para compatibilidad con declustering
if 'magnitude' in df_base.columns and 'mag' not in df_base.columns:
    df_base = df_base.rename(columns={'magnitude': 'mag'})

print(f'Catalogo base: N = {len(df_base)} eventos')
print(f'Rango temporal: {df_base.time.min().date()} a {df_base.time.max().date()}')
print(f'Magnitud: min = {df_base.mag.min():.1f}, max = {df_base.mag.max():.1f}')
df_base[['time', 'latitude', 'longitude', 'depth', 'mag', 'dip1', 'dip2']].head()\
"""))

cells.append(nbf.v4.new_markdown_cell("""\
## 2. Aplicacion de los tres metodos de declustering

- **Ninguno**: catalogo intacto (control negativo: incluye replicas).
- **GK fijo**: ventanas (100 km, 30 dias). Replica del enfoque de NB10.
- **GK variable**: ventanas R(M), T(M) crecientes con magnitud (literatura estandar).
"""))

cells.append(nbf.v4.new_code_cell("""\
catalogos = {
    'none':        df_base.copy(),
    'gk_fixed':    gardner_knopoff_fixed(df_base, radius_km=100.0, window_days=30.0),
    'gk_variable': gardner_knopoff_variable(df_base),
}

print('| Metodo         | N   | Reduccion |')
print('|----------------|-----|-----------|')
for name, df in catalogos.items():
    pct = 100 * (1 - len(df) / len(df_base))
    print(f'| {name:<14} | {len(df):>3} | {pct:>6.1f}%  |')\
"""))

cells.append(nbf.v4.new_markdown_cell("""\
## 3. Monte Carlo 1000x para cada catalogo

`time_shuffling_null` pre-computa la matriz NxN de Coulomb stress con UNA sola
llamada vectorizada a Skyfield, luego ejecuta 1000 permutaciones en paralelo.

Esto es ~10x mas rapido que el loop original de NB10.
"""))

cells.append(nbf.v4.new_code_cell("""\
N_ITER = 1000
MU = 0.4

resultados = {}
for name, df in catalogos.items():
    t0 = time.perf_counter()
    res = time_shuffling_null(df, n_iterations=N_ITER, mu=MU, seed=42, n_jobs=-1)
    res['elapsed_s'] = time.perf_counter() - t0
    resultados[name] = res
    print(f'[{name:<14}] N={res[\"n_events\"]:>3}  '
          f'obs={res[\"observed_fraction\"]:.4f}  '
          f'null={res[\"null_mean\"]:.4f}+-{res[\"null_std\"]:.4f}  '
          f'p={res[\"p_value\"]:.4f}  '
          f'({res[\"elapsed_s\"]:.1f}s)')\
"""))

cells.append(nbf.v4.new_markdown_cell("""\
## 4. Matriz de sensibilidad final

Tabla compacta para discusion en el manuscrito: si los tres metodos producen
p-values en el mismo orden de magnitud, la conclusion es ROBUSTA. Si solo el
catalogo sin declusterizar da p<0.05 mientras GK_variable da p>0.5, la senal
era artefacto de las replicas (= triggering de Coulomb estatico, no marea).
"""))

cells.append(nbf.v4.new_code_cell("""\
tabla = pd.DataFrame({
    name: {
        'N': res['n_events'],
        'observado_pct':   100 * res['observed_fraction'],
        'null_mean_pct':   100 * res['null_mean'],
        'null_std_pct':    100 * res['null_std'],
        'exceso_pct':      100 * (res['observed_fraction'] - res['null_mean']),
        'p_value':         res['p_value'],
    }
    for name, res in resultados.items()
}).T

# Guarda tabla para el manuscrito
tabla_path = RESULTS_STATS / 'declustering_sensitivity.csv'
tabla.to_csv(tabla_path)
print(f'Guardado: {tabla_path.relative_to(tabla_path.parents[2])}')
tabla\
"""))

cells.append(nbf.v4.new_markdown_cell("""\
## 5. Visualizacion: distribuciones nulas vs observado
"""))

cells.append(nbf.v4.new_code_cell("""\
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
colors = {'none': '#d62728', 'gk_fixed': '#ff7f0e', 'gk_variable': '#2ca02c'}

for ax, (name, res) in zip(axes, resultados.items()):
    null = res['null_distribution']
    obs = res['observed_fraction']
    ax.hist(null * 100, bins=40, color=colors[name], alpha=0.7, edgecolor='black')
    ax.axvline(obs * 100, color='black', linestyle='--', linewidth=2,
               label=f'Observado = {obs*100:.2f}%')
    ax.set_title(f'{name}  (N={res[\"n_events\"]}, p={res[\"p_value\"]:.3f})')
    ax.set_xlabel('Fraccion CFS > 0  (%)')
    ax.legend(loc='upper right', fontsize=9)

axes[0].set_ylabel('Frecuencia (de 1000 iter)')
fig.suptitle('Distribucion nula MC bajo distintos metodos de declustering', fontsize=13)
fig.tight_layout()

fig_path = RESULTS_FIGURES / 'declustering_null_distributions.png'
fig.savefig(fig_path, dpi=150, bbox_inches='tight')
print(f'Guardado: {fig_path.relative_to(fig_path.parents[2])}')
plt.show()\
"""))

cells.append(nbf.v4.new_markdown_cell("""\
## 6. Conclusion

La robustez del *Lunar Tidal Triggering* depende del comportamiento del p-value
bajo declustering progresivamente mas estricto:

- Si **p crece monotonicamente** al declusterizar (`none` -> `gk_variable`),
  la senal original incluia agrupamientos espurios por replicas y no es real.
- Si **p se mantiene estable o solo crece moderadamente**, la senal es robusta
  y refleja modulacion gravitacional de eventos independientes.

Los valores numericos exactos quedan en `results/statistics/declustering_sensitivity.csv`
para citas y revision por pares.
"""))

nb.cells = cells
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.10"},
}

out_path = Path("notebooks") / "11_declustering_sensitivity.ipynb"
nbf.write(nb, out_path)
print(f'Escrito: {out_path}')
