import nbformat as nbf
import os

# Create the notebook object
nb = nbf.v4.new_notebook()

cells = []

# Cell 1: Markdown Title
cells.append(nbf.v4.new_markdown_cell("""# Notebook 05: Análisis de Retraso Elástico (Lag) y Derivadas
**Autor:** Iván Andrés Mena Contreras
**Proyecto:** Lunar Tidal Triggering of Earthquakes

En el cuaderno anterior concluimos que el gatillo instantáneo (la fase lunar en t=0) es indistinguible del ruido ($p=0.313$). 

**Objetivo:** Investigar si el gatillo lunar opera con un retraso temporal debido a la viscoelasticidad de la corteza (retraso en la transferencia del estrés de marea a la falla tectónica), o si el factor detonante depende de la tasa de cambio orbital de la marea ($dF/dt$), es decir, si los sismos se agrupan cuando la fuerza gravitacional varía más drásticamente."""))

# Cell 2: Imports and Data Loading
cells.append(nbf.v4.new_code_cell("""import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astroquery.jplhorizons import Horizons
from astropy.time import Time
from tqdm import tqdm
import warnings
import os

warnings.filterwarnings('ignore')

# Configuración visual
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
os.makedirs('../results/figures', exist_ok=True)

# Carga de datos del subset ecuatorial
data_path = '../data/processed/earthquakes_ecuatorial_subset.csv'
df = pd.read_csv(data_path)
df['time'] = pd.to_datetime(df['time'], format='mixed', utc=True)

print(f"Catálogo cargado: {len(df)} eventos del cinturón ecuatorial superficial.")"""))

# Cell 3: Cálculo de dF/dt
cells.append(nbf.v4.new_markdown_cell("""## 1. Tasa de Cambio Orbital ($dF/dt$)
Como la órbita de la Luna es elíptica, la fuerza de marea no solo varía por la fase, sino también por el perigeo y apogeo. Para capturar este componente dinámico, necesitamos saber si la Luna se acercaba o se alejaba y a qué velocidad (Aproximación de derivada: $\\Delta \\text{Distancia} / \\Delta \\text{Tiempo}$).

Consultamos la distancia a Horizons 2 horas antes de cada evento y calculamos la tasa de cambio."""))

cells.append(nbf.v4.new_code_cell("""# Inicializamos columnas
df['dist_minus_2h'] = np.nan

print("Consultando JPL Horizons para t - 2h...")
for idx, row in tqdm(df.iterrows(), total=len(df)):
    try:
        t_0 = Time(row['time'])
        t_minus2 = t_0 - 2/24.0 # 2 horas antes (en días)
        
        # Consultamos la luna (id=301) desde el centro de la Tierra (500)
        obj = Horizons(id='301', location='500', epochs=t_minus2.jd)
        eph = obj.ephemerides()
        df.loc[idx, 'dist_minus_2h'] = eph['delta'][0] # distancia en AU
    except Exception as e:
        print(f"Error en índice {idx}: {e}")

# Tasa de cambio: delta(AU) / 2h.
# Unidades: (AU / hr) * 10^6 para legibilidad.
df['dist_rate_of_change'] = ((df['moon_distance_au'] - df['dist_minus_2h']) / 2.0) * 1e6

print("Cálculo de dF/dt completado.")
df[['time', 'moon_distance_au', 'dist_minus_2h', 'dist_rate_of_change']].head()"""))

# Cell 4: Histograma de Derivadas
cells.append(nbf.v4.new_markdown_cell("""## 2. Distribución de la Tasa de Cambio
Visualizamos la distribución de `dist_rate_of_change`. 
- Valores negativos: La Luna se estaba acercando a la Tierra.
- Valores positivos: La Luna se estaba alejando de la Tierra.
- Valores cercanos a cero: La Luna estaba estacionaria en Perigeo o Apogeo."""))

cells.append(nbf.v4.new_code_cell("""plt.figure(figsize=(10, 6))
sns.histplot(df['dist_rate_of_change'].dropna(), bins=15, kde=True, color='#2ecc71')
plt.axvline(0, color='red', linestyle='--', linewidth=2, label='Estacionario (Perigeo/Apogeo)')
plt.title('Distribución de la Tasa de Cambio de la Distancia Lunar ($dF/dt$)', fontsize=14)
plt.xlabel('Tasa de Cambio ($\\mu$AU / hora)')
plt.ylabel('Frecuencia')
plt.legend()
plt.tight_layout()

dist_path = '../results/figures/05_derivative_dist.png'
plt.savefig(dist_path, dpi=300)
print(f"Histograma guardado en: {dist_path}")
plt.show()"""))

# Cell 5: Lag Analysis (Fases)
cells.append(nbf.v4.new_markdown_cell("""## 3. Retraso Viscoelástico (Lag Analysis)
Extraemos el ángulo de fase lunar (RA) exactamente 6 horas antes ($t-6h$) y 12 horas antes ($t-12h$) para cada evento, usando JPL Horizons."""))

cells.append(nbf.v4.new_code_cell("""df['moon_ra_minus_6h'] = np.nan
df['moon_ra_minus_12h'] = np.nan

print("Consultando JPL Horizons para t-6h y t-12h...")
for idx, row in tqdm(df.iterrows(), total=len(df)):
    try:
        t_0 = Time(row['time'])
        t_minus6 = t_0 - 6/24.0
        t_minus12 = t_0 - 12/24.0
        
        # Consulta en batch (lista de épocas)
        obj = Horizons(id='301', location='500', epochs=[t_minus6.jd, t_minus12.jd])
        eph = obj.ephemerides()
        
        df.loc[idx, 'moon_ra_minus_6h'] = eph['RA'][0]
        df.loc[idx, 'moon_ra_minus_12h'] = eph['RA'][1]
    except Exception as e:
        print(f"Error en índice {idx}: {e}")

# Ajuste temporal simple para el test de fase (Fase = RA Luna - RA Sol)
# Simplificamos asumiendo que el desplazamiento del sol en 6h es mínimo, 
# la diferencia principal radica en el vector RA de la luna en relación con la rotación.
# Para el test rotacional direccional usamos directamente moon_ra de manera comparable.

print("Extracción temporal completada.")"""))

# Cell 6: Schuster Lag Tests
cells.append(nbf.v4.new_markdown_cell("""## 4. Test de Schuster con Lags
Ejecutamos el Test de Schuster para evaluar si la señal se agrupa significativamente en el pasado."""))

cells.append(nbf.v4.new_code_cell("""def schuster_test(angulos_deg):
    \"\"\"Calcula el p-value del Test de Schuster basado en ángulos en grados.\"\"\"
    angulos_deg = angulos_deg.dropna()
    angulos_rad = np.radians(angulos_deg)
    N = len(angulos_rad)
    if N == 0: return np.nan, np.nan
    
    R = np.sqrt(np.sum(np.cos(angulos_rad))**2 + np.sum(np.sin(angulos_rad))**2)
    p_value = np.exp(-(R**2) / N)
    return R, p_value

# La hipótesis original usa la fase, que se mapea directamente sobre moon_ra
# como el vector principal de influencia ecuatorial
R_0, p_0 = schuster_test(df['moon_ra'])
R_6, p_6 = schuster_test(df['moon_ra_minus_6h'])
R_12, p_12 = schuster_test(df['moon_ra_minus_12h'])

print("=== TEST DE SCHUSTER: ANÁLISIS DE RETRASO VISCOELÁSTICO ===")
print(f"t = 0h  (Instantáneo): p-value = {p_0:.5f}")
print(f"t = -6h (Retraso 6h):  p-value = {p_6:.5f}")
print(f"t = -12h(Retraso 12h): p-value = {p_12:.5f}")

if min(p_0, p_6, p_12) < 0.05:
    print("\\n-> ¡SEÑAL DETECTADA! La corteza responde con retardo y muestra agrupamiento estadístico.")
else:
    print("\\n-> CONCLUSIÓN FINAL: Ni el estrés instantáneo ni el retrasado viscoelástico logran rechazar la hipótesis nula ($p > 0.05$).")"""))

# Cell 7: Export
cells.append(nbf.v4.new_markdown_cell("""## 5. Exportar Dataset
Guardamos la matriz final con todos los derivados computados."""))

cells.append(nbf.v4.new_code_cell("""export_path = '../data/processed/earthquakes_lag_analysis.csv'
df.to_csv(export_path, index=False)
print(f"Dataset final de lags exportado a: {export_path}")"""))

# Assign cells to notebook
nb['cells'] = cells

# Save the notebook
notebook_path = r'c:\Users\IVAN MENA\Documents\lunar-seismic-triggering\notebooks\05_lag_analysis.ipynb'
with open(notebook_path, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"Notebook created at {notebook_path}")
