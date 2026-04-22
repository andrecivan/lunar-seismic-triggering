import nbformat as nbf
import os

# Create the notebook object
nb = nbf.v4.new_notebook()

cells = []

# Cell 1: Markdown Title
cells.append(nbf.v4.new_markdown_cell("""# Notebook 04: Visualización Espacial y Segmentación Ecuatorial
**Autor:** Iván Andrés Mena Contreras
**Proyecto:** Lunar Tidal Triggering of Earthquakes

En el cuaderno anterior confirmamos que, a escala global, los grandes sismos ($M \geq 7.0$) se distribuyen como ruido respecto a las fases lunares ($p \approx 0.45$). 

**Objetivo de este Cuaderno:**
Aislar la señal lunar cruzando la fase lunar con la ubicación geográfica. Al calcular el Punto Sublunar y segmentar geográficamente (específicamente en el cinturón de subducción ecuatorial, entre latitudes -30° y +30°), buscamos detectar si la marea de la corteza, que es máxima en el ecuador terrestre, actúa como gatillo regional."""))

# Cell 2: Imports and Data Loading
cells.append(nbf.v4.new_code_cell("""import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import folium
import os
from astropy.time import Time
from astropy.coordinates import get_sun, get_body, EarthLocation, AltAz
from astropy import units as u

# Configuración visual
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)

# Crear directorios
os.makedirs('../results/maps', exist_ok=True)
os.makedirs('../results/figures', exist_ok=True)

# Carga de datos limpios
data_path = '../data/processed/earthquakes_with_moon_final.csv'
df = pd.read_csv(data_path)
df['time'] = pd.to_datetime(df['time'], format='mixed')

print(f"Catálogo cargado: {len(df)} eventos.")"""))

# Cell 3: Global Map
cells.append(nbf.v4.new_markdown_cell("""## 1. Visualización Espacial Global
Graficamos todos los epicentros para comprender la distribución tectónica de nuestro catálogo."""))

cells.append(nbf.v4.new_code_cell("""# Función para determinar color según profundidad
def color_profundidad(prof):
    if prof <= 20: return '#e74c3c'  # Rojo (Superficial)
    elif prof <= 45: return '#f39c12' # Naranja (Intermedio)
    return '#2980b9'                  # Azul (Profundo)

# Crear mapa base
mapa = folium.Map(location=[0, 0], zoom_start=2, tiles='CartoDB dark_matter')

for idx, row in df.iterrows():
    # El radio depende de la magnitud
    radio = np.exp(row['magnitude'] - 6.5) * 2
    
    folium.CircleMarker(
        location=[row['latitude'], row['longitude']],
        radius=radio,
        color=color_profundidad(row['depth']),
        fill=True,
        fill_color=color_profundidad(row['depth']),
        fill_opacity=0.7,
        tooltip=f"M{row['magnitude']} - {row['place']} (Prof: {row['depth']}km)"
    ).add_to(mapa)

map_path = '../results/maps/04_epicenters_map.html'
mapa.save(map_path)
print(f"Mapa global guardado en: {map_path}")
mapa"""))

# Cell 4: Sublunar Point
cells.append(nbf.v4.new_markdown_cell("""## 2. Cálculo del Punto Sublunar
El "Punto Sublunar" es el lugar geográfico exacto en la Tierra donde la Luna se encontraba exactamente en el cenit (arriba en el cielo) en el momento preciso de la ruptura del sismo.

La Latitud del punto sublunar es igual a la Declinación de la Luna (`moon_dec`).
La Longitud depende del Tiempo Sideral Aparente de Greenwich (GAST) y la Ascensión Recta (`moon_ra`)."""))

cells.append(nbf.v4.new_code_cell("""# Astropy Time object
tiempos = Time(df['time'])

# Calcular el GAST (Greenwich Apparent Sidereal Time)
gast = tiempos.sidereal_time('apparent', 'greenwich')

# La longitud del punto sublunar es: moon_ra - gast
# (Convirtiendo todo a grados y normalizando entre -180 y +180)
lon_sublunar = (df['moon_ra'] - gast.degree + 180) % 360 - 180

# Asignar al DataFrame
df['sublunar_lat'] = df['moon_dec']
df['sublunar_lon'] = lon_sublunar

print("Cálculo del punto sublunar completado.")
df[['time', 'latitude', 'longitude', 'sublunar_lat', 'sublunar_lon']].head()"""))

# Cell 5: Angular Distance
cells.append(nbf.v4.new_markdown_cell("""## 3. Distancia Angular al Punto Sublunar
Calculamos la distancia ortodrómica (o ángulo central) entre el epicentro y el punto sublunar usando trigonometría esférica."""))

cells.append(nbf.v4.new_code_cell("""def distancia_angular(lat1, lon1, lat2, lon2):
    \"\"\"Calcula la distancia angular en grados entre dos puntos en una esfera.\"\"\"
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    
    # Ley esférica de los cosenos
    d_sigma = np.arccos(
        np.sin(lat1) * np.sin(lat2) + 
        np.cos(lat1) * np.cos(lat2) * np.cos(abs(lon1 - lon2))
    )
    return np.degrees(d_sigma)

# Calcular distancias
df['angular_distance'] = distancia_angular(
    df['latitude'], df['longitude'], 
    df['sublunar_lat'], df['sublunar_lon']
)

# Visualizar distribución
plt.figure(figsize=(10, 6))
sns.histplot(df['angular_distance'], bins=20, kde=True, color='#8e44ad')
plt.axvline(90, color='red', linestyle='--', label='Plano del Horizonte (90°)')
plt.title('Distribución de la Distancia Angular entre Epicentro y Punto Sublunar', fontsize=14)
plt.xlabel('Distancia Angular (Grados)')
plt.ylabel('Frecuencia')
plt.legend()
plt.tight_layout()

dist_path = '../results/figures/04_angular_dist.png'
plt.savefig(dist_path, dpi=300)
print(f"Histograma guardado en: {dist_path}")
plt.show()"""))

# Cell 6: Geographic Filter
cells.append(nbf.v4.new_markdown_cell("""## 4. Filtro Ecuatorial Regional
La fuerza de marea sólida de la Tierra es mayor en el ecuador. Vamos a aislar nuestra búsqueda de señal filtrando exclusivamente sismos superficiales ($<20$ km) dentro de las latitudes bajas ($-30^{\circ}$ a $+30^{\circ}$)."""))

cells.append(nbf.v4.new_code_cell("""# Filtrar Sismos Superficiales y Ecuatoriales
condicion_profundidad = df['depth'] < 20
condicion_latitud = (df['latitude'] >= -30) & (df['latitude'] <= 30)

df_ecuatorial = df[condicion_profundidad & condicion_latitud].copy()

print(f"Sismos originales en el catálogo: {len(df)}")
print(f"Sismos en el Cinturón Ecuatorial Superficial: {len(df_ecuatorial)}")"""))

# Cell 7: Regional Schuster Test
cells.append(nbf.v4.new_markdown_cell("""## 5. Test de Schuster Regional
Recalculamos el valor $p$ utilizando únicamente este nuevo sub-grupo altamente susceptible a la fuerza de marea cortical máxima."""))

cells.append(nbf.v4.new_code_cell("""def schuster_test(angulos_deg):
    \"\"\"Calcula el p-value del Test de Schuster basado en ángulos en grados.\"\"\"
    angulos_rad = np.radians(angulos_deg)
    N = len(angulos_rad)
    if N == 0: return np.nan, np.nan
    
    R = np.sqrt(np.sum(np.cos(angulos_rad))**2 + np.sum(np.sin(angulos_rad))**2)
    p_value = np.exp(-(R**2) / N)
    return R, p_value

# Extraer el ángulo de fase lunar que usamos en el Nb 03
# Fase = (RA Luna - RA Sol)
sol = get_sun(Time(df_ecuatorial['time']))
luna = get_body('moon', Time(df_ecuatorial['time']))
fase_lunar = np.mod(luna.ra.degree - sol.ra.degree, 360)

# Calcular Schuster
R_reg, p_reg = schuster_test(fase_lunar)

print("=== TEST DE SCHUSTER REGIONAL (ECUATORIAL SUPERFICIAL) ===")
print(f"  N = {len(df_ecuatorial)}")
print(f"  Vector R: {R_reg:.4f}")
print(f"  p-value:  {p_reg:.5f}")

if p_reg < 0.05:
    print("  -> ¡CONCLUSIÓN: Se rechaza la hipótesis nula! Hay una correlación ESTADÍSTICAMENTE SIGNIFICATIVA con la fase lunar en esta región.")
else:
    print("  -> CONCLUSIÓN: No se rechaza la hipótesis nula. El subgrupo regional sigue comportándose como ruido estadístico.")"""))

# Cell 8: Export
cells.append(nbf.v4.new_markdown_cell("""## 6. Exportar Datos
Guardamos este subset de alto interés para posibles cruces de datos con fallas tectónicas específicas en el siguiente cuaderno."""))

cells.append(nbf.v4.new_code_cell("""export_path = '../data/processed/earthquakes_ecuatorial_subset.csv'
df_ecuatorial.to_csv(export_path, index=False)
print(f"Subset exportado exitosamente a: {export_path}")"""))

# Assign cells to notebook
nb['cells'] = cells

# Save the notebook
notebook_path = r'c:\Users\IVAN MENA\Documents\lunar-seismic-triggering\notebooks\04_spatial_analysis.ipynb'
with open(notebook_path, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"Notebook created at {notebook_path}")
