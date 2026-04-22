import nbformat as nbf
import os

# Create the notebook object
nb = nbf.v4.new_notebook()

cells = []

# Cell 1: Markdown Title
cells.append(nbf.v4.new_markdown_cell("""# Cuaderno 03: Análisis Estadístico Exploratorio y Test de Schuster
**Autor:** Iván Andrés Mena Contreras
**Proyecto:** Lunar Tidal Triggering of Earthquakes

En este cuaderno realizaremos el análisis estadístico del catálogo filtrado ($M \geq 7.0$). El objetivo es identificar si existe una correlación estadística significativa entre las fases lunares y la ocurrencia de sismos, utilizando estadística circular y el Test de Schuster.

## Objetivos
1. **Cálculo de la Fase Lunar:** Estimar el ángulo de fase sinódica lunar (0°-360°) para cada sismo.
2. **Análisis Circular y Visualización:** Mapear la distribución temporal de eventos en un plano polar según la profundidad.
3. **Test de Schuster:** Evaluar estadísticamente si los sismos ocurren con una distribución uniforme (hipótesis nula) o si hay agrupación preferencial."""))

# Cell 2: Imports
cells.append(nbf.v4.new_code_cell("""import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.time import Time
from astropy.coordinates import get_sun, get_body
from astropy import units as u

# Configuración de estilos visuales
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 12

# Crear directorios de salida si no existen
os.makedirs('../results/figures', exist_ok=True)"""))

# Cell 3: Carga de Datos
cells.append(nbf.v4.new_markdown_cell("""## 1. Carga y Preparación de Datos
Cargamos el catálogo previamente procesado."""))

cells.append(nbf.v4.new_code_cell("""# Cargar datos
data_path = '../data/processed/earthquakes_with_moon_final.csv'
df = pd.read_csv(data_path)

# Convertir time a datetime
df['time'] = pd.to_datetime(df['time'], format='ISO8601')

print(f"Catálogo cargado: {len(df)} eventos")
df.head()"""))

# Cell 4: Cálculo de Fase Lunar
cells.append(nbf.v4.new_markdown_cell("""## 2. Cálculo de la Fase Lunar
Para calcular la fase lunar de manera precisa, obtenemos las coordenadas ecuatoriales del Sol y la Luna utilizando Astropy. La fase puede aproximarse como la diferencia en Ascensión Recta (RA) entre la Luna y el Sol.

Una fase de $0^{\circ}$ corresponde a Luna Nueva, $90^{\circ}$ a Cuarto Creciente, $180^{\circ}$ a Luna Llena y $270^{\circ}$ a Cuarto Menguante."""))

cells.append(nbf.v4.new_code_cell("""def calcular_fase_lunar(fechas):
    \"\"\"
    Calcula la fase lunar en grados (0 a 360) para una serie de fechas.
    \"\"\"
    tiempos = Time(fechas)
    
    # Obtener posiciones aparentes del Sol y la Luna
    sol = get_sun(tiempos)
    luna = get_body('moon', tiempos)
    
    # La elongación (diferencia en Ascensión Recta)
    # Se ajusta para que el rango sea de 0 a 360 grados
    elongacion = luna.ra.degree - sol.ra.degree
    fase = np.mod(elongacion, 360.0)
    
    return fase

# Aplicar el cálculo al dataset
print("Calculando fases lunares...")
df['moon_phase_deg'] = calcular_fase_lunar(df['time'])
df['moon_phase_rad'] = np.radians(df['moon_phase_deg'])

print("Fases calculadas. Estadísticas:")
print(df['moon_phase_deg'].describe())"""))

# Cell 5: Segmentación por Profundidad
cells.append(nbf.v4.new_markdown_cell("""## 3. Segmentación por Profundidad
Se han observado diferencias en los mecanismos de ruptura de acuerdo a la profundidad. Separamos el catálogo en tres rangos (bins) corticales para análisis diferenciado:
- Superficiales: $0 - 20$ km
- Intermedios: $20 - 45$ km
- Profundos: $45 - 70$ km (y más)"""))

cells.append(nbf.v4.new_code_cell("""def categorizar_profundidad(depth):
    if depth <= 20:
        return 'Superficial (0-20 km)'
    elif depth <= 45:
        return 'Intermedio (20-45 km)'
    else:
        return 'Profundo (>45 km)'

df['depth_category'] = df['depth'].apply(categorizar_profundidad)

print("Distribución por profundidad:")
print(df['depth_category'].value_counts())"""))

# Cell 6: Visualización Polar
cells.append(nbf.v4.new_markdown_cell("""## 4. Visualización Polar de las Fases
Visualizamos la distribución de sismos en un gráfico polar para identificar patrones de concentración. Cada punto representa un evento, posicionado según su fase lunar, tamaño según su magnitud y color según su profundidad."""))

cells.append(nbf.v4.new_code_cell("""fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})

colores = {
    'Superficial (0-20 km)': '#e74c3c',  # Rojo
    'Intermedio (20-45 km)': '#f39c12',  # Naranja
    'Profundo (>45 km)': '#3498db'       # Azul
}

for categoria, color in colores.items():
    subset = df[df['depth_category'] == categoria]
    
    # El tamaño del marcador escalado exponencialmente según la magnitud
    tamanos = np.exp(subset['magnitude'] - 6.5) * 20
    
    # Ploteamos en el eje polar (theta=radianes, r=magnitud)
    # R se fijará alrededor de un radio base para no dispersar demasiado por magnitud,
    # pero podemos hacer que r refleje una constante o la magnitud misma.
    # Usaremos magnitud como radio para mayor claridad visual.
    ax.scatter(subset['moon_phase_rad'], subset['magnitude'], 
               c=color, s=tamanos, alpha=0.7, edgecolors='k', label=categoria)

# Configuración de los ejes polares
ax.set_theta_zero_location('N') # 0° arriba (Luna Nueva)
ax.set_theta_direction(-1) # Sentido horario
ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
ax.set_xticklabels(['Nueva (0°)', '45°', 'Creciente (90°)', '135°', 
                    'Llena (180°)', '225°', 'Menguante (270°)', '315°'])

# Ajuste del radio
ax.set_ylim(6.5, 9.5)
ax.set_yticks([7.0, 8.0, 9.0])
ax.set_yticklabels(['M7', 'M8', 'M9'], color='grey')
ax.set_title('Distribución Polar de Terremotos (M$\\geq$7.0) por Fase Lunar', fontsize=16, pad=20)
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

plt.tight_layout()
plt.savefig('../results/figures/03_polar_distribution.png', dpi=300, bbox_inches='tight')
plt.show()"""))

# Cell 7: Test de Schuster
cells.append(nbf.v4.new_markdown_cell("""## 5. Test de Schuster
El Test de Schuster es una prueba estadística rigurosa para evaluar la periodicidad en series temporales direccionales.
Calculamos el vector suma $R$:
$R = \sqrt{(\sum \cos \phi_i)^2 + (\sum \sin \phi_i)^2}$

La probabilidad $p$-value de que esta distribución provenga de una distribución aleatoria uniforme se aproxima como:
$p = e^{-R^2/N}$

Si $p < 0.05$, rechazamos la hipótesis nula y concluimos que existe una periodicidad o agrupamiento significativo."""))

cells.append(nbf.v4.new_code_cell("""def schuster_test(fases_rad):
    \"\"\"
    Realiza el Test de Schuster para detectar agrupamientos circulares.
    Retorna el vector R y el valor p.
    \"\"\"
    N = len(fases_rad)
    if N == 0:
        return np.nan, np.nan
        
    sum_cos = np.sum(np.cos(fases_rad))
    sum_sin = np.sum(np.sin(fases_rad))
    
    R = np.sqrt(sum_cos**2 + sum_sin**2)
    p_value = np.exp(-(R**2) / N)
    
    return R, p_value

# 1. Test para todo el catálogo
R_total, p_total = schuster_test(df['moon_phase_rad'])

# 2. Test para sismos superficiales (<20km)
sismos_superficiales = df[df['depth'] < 20]
R_superf, p_superf = schuster_test(sismos_superficiales['moon_phase_rad'])

print("=== RESULTADOS DEL TEST DE SCHUSTER ===")
print(f"Catálogo Completo (N={len(df)}):")
print(f"  Vector R: {R_total:.4f}")
print(f"  p-value:  {p_total:.4g} " + ("(Significativo)" if p_total < 0.05 else "(No significativo)"))
print()
print(f"Sismos Superficiales <20km (N={len(sismos_superficiales)}):")
print(f"  Vector R: {R_superf:.4f}")
print(f"  p-value:  {p_superf:.4g} " + ("(Significativo)" if p_superf < 0.05 else "(No significativo)"))"""))

# Cell 8: Histograma y Distribución
cells.append(nbf.v4.new_markdown_cell("""## 6. Histograma de Frecuencias
Finalmente, observaremos la distribución en un histograma con bins de $15^{\circ}$ para verificar visualmente cualquier desviación respecto a una distribución uniforme (hipótesis nula)."""))

cells.append(nbf.v4.new_code_cell("""fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Bins de 15 grados
bins = np.arange(0, 361, 15)

# Plot para catálogo completo
counts_total, _, _ = ax1.hist(df['moon_phase_deg'], bins=bins, color='#7f8c8d', alpha=0.7, edgecolor='black')
ax1.axhline(len(df) / (360/15), color='red', linestyle='--', label='Distribución Uniforme Esperada')
ax1.set_title('Distribución de Frecuencias: Catálogo Completo (M$\\geq$7.0)', fontsize=14)
ax1.set_ylabel('Número de Eventos')
ax1.legend()
ax1.grid(axis='y', alpha=0.3)

# Plot para superficiales
counts_sup, _, _ = ax2.hist(sismos_superficiales['moon_phase_deg'], bins=bins, color='#e74c3c', alpha=0.7, edgecolor='black')
ax2.axhline(len(sismos_superficiales) / (360/15), color='red', linestyle='--', label='Distribución Uniforme Esperada')
ax2.set_title('Distribución de Frecuencias: Sismos Superficiales (<20km)', fontsize=14)
ax2.set_xlabel('Fase Lunar (grados)')
ax2.set_ylabel('Número de Eventos')
ax2.set_xticks(np.arange(0, 361, 45))
ax2.set_xticklabels(['0°\\n(Nueva)', '45°', '90°\\n(Creciente)', '135°', '180°\\n(Llena)', '225°', '270°\\n(Menguante)', '315°', '360°'])
ax2.legend()
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('../results/figures/03_schuster_test.png', dpi=300, bbox_inches='tight')
plt.show()"""))

# Assign cells to notebook
nb['cells'] = cells

# Save the notebook
notebook_path = r'c:\Users\IVAN MENA\Documents\lunar-seismic-triggering\notebooks\03_statistical_analysis.ipynb'
with open(notebook_path, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"Notebook created at {notebook_path}")
