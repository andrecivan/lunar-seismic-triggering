# Conclusiones Finales — Revisiting Lunar Tidal Triggering

**Autor:** Ing. Iván Andrés Mena Contreras
**Fecha de revisión:** 2026-04-21
**Estado:** Hipótesis rechazada tras análisis de robustez. *No se puede rechazar la hipótesis nula.*

---

## 1. Punto de quiebre: el p-value bajo declustering progresivo

La señal estadística reportada en versiones preliminares del manuscrito (Notebook 09, sin declustering: $p = 0.025$) **no sobrevive** al someter el catálogo a procedimientos estándar de eliminación de réplicas. La tabla siguiente documenta el deterioro monotónico del soporte estadístico a medida que el declustering se vuelve más estricto:

| Procedimiento | Notebook | N (eventos) | % $\Delta$CFS > 0 (real) | Media nula MC | **p-value empírico** | Decisión a $\alpha = 0.05$ |
|---|---|---|---|---|---|---|
| **Sin declustering** | 09 | 333 | 64.86 % | 60.39 % | **0.025** | Rechaza H₀ |
| **Gardner-Knopoff fijo** (100 km / 30 d) | 10 | 309 | 63.75 % | 59.65 % | **0.050** | Borderline (umbral exacto) |
| **Gardner-Knopoff variable** (Helmstetter & Sornette 2003) | 11 | 288 | 63.19 % | 59.42 % | **0.080** | **No rechaza H₀** |

Lectura: cada paso de declustering (a) reduce el tamaño efectivo del catálogo, (b) reduce el exceso observado de eventos favorecidos por la marea sobre la línea base nula, y (c) **eleva el p-value empírico**.

### 1.1. Localización del cruce
El umbral de significancia $\alpha = 0.05$ es cruzado **exactamente al aplicar Gardner-Knopoff fijo** (NB10: $p = 0.0500$, frontera exacta). Bajo el procedimiento parametrizado por magnitud (Helmstetter & Sornette), la señal se aleja sin ambigüedad de la región de rechazo ($p = 0.080$, equivalente a una desviación de menos de $1.5\sigma$ respecto al nulo).

### 1.2. Implicación
La aparente significancia del análisis sin declustering es **atribuible al agrupamiento espacio-temporal de réplicas**, no a una influencia gravitacional lunar genuina. Las réplicas comparten una causa física común (transferencia de Coulomb estática del mainshock asociado) y su presencia en el catálogo viola la hipótesis de independencia temporal sobre la que se construye el modelo nulo de time-shuffling.

---

## 2. El sesgo geométrico del 60 %

Una observación adicional, producto colateral de la simulación Monte Carlo, merece atención metodológica explícita: **la línea base bajo la hipótesis nula no se centra en el 50 %**, como un razonamiento ingenuo sobre simetría sugeriría. En las tres configuraciones del catálogo (333, 309 y 288 eventos), la fracción media de eventos con $\Delta$CFS > 0 bajo time-shuffling oscila en torno al **60 %**:

$$
\bar{p}_{nulo} \approx 60\% \quad (\text{rango: } 59.4\% - 60.4\%)
$$

### 2.1. Origen del sesgo

Este desplazamiento de la línea base no es un artefacto numérico ni un error de implementación; es una propiedad **geométrica intrínseca** del catálogo de mega-terremotos terrestres. Surge de la conjunción de tres asimetrías:

1. **Distribución latitudinal de las fallas:** las zonas de subducción que producen $M \geq 7$ se concentran en el Cinturón de Fuego del Pacífico (latitudes medias y bajas, predominantemente en el hemisferio norte).
2. **Geometría de subducción:** los planos de falla activos tienen ángulos de buzamiento ($dip$) acotados (típicamente $20° - 60°$); no son una distribución uniforme en $[0°, 90°]$.
3. **Declinación lunar:** la órbita lunar mantiene a la Luna oscilando aproximadamente entre $\pm 28.5°$ de declinación; sus posiciones aparentes desde las zonas sísmicas privilegian configuraciones geométricas particulares.

La interacción de estas tres asimetrías produce que **incluso bajo asignación temporal aleatoria**, la geometría preferencial del sistema Tierra-Luna favorece configuraciones de $\Delta$CFS > 0 con probabilidad cercana a 0.60, no a 0.50.

### 2.2. Consecuencia metodológica

Cualquier estudio futuro de triggering tidal **debe** estimar empíricamente la línea base nula mediante simulación Monte Carlo específica del catálogo y la geometría usada. **Comparar contra el 50 % es estadísticamente incorrecto** y produce conclusiones falsamente positivas: el "exceso" del 64.86 % observado en NB09 es solo $\approx 4.5$ puntos porcentuales sobre la línea geométrica del 60.39 %, no $\approx 14.9$ puntos sobre el 50 % ingenuo.

Este hallazgo es, paradójicamente, el resultado científicamente más sólido del estudio: una **caracterización cuantitativa del sesgo geométrico de fondo** que aparece al combinar el catálogo USGS de mega-sismos con la órbita lunar.

---

## 3. Veredicto final

> **No se puede rechazar la hipótesis nula** ($H_0$: la marea gravitacional lunar no modula la ocurrencia de mega-terremotos $M \geq 7$).

### 3.1. Justificación
Los datos disponibles, una vez sometidos a:
- Resolución desesgada de la ambigüedad nodal del mecanismo focal (promedio de los dos planos),
- Eliminación de réplicas mediante el método estándar de Gardner-Knopoff con ventanas dependientes de magnitud,
- Calibración empírica de la línea base geométrica vía Monte Carlo time-shuffling con 1000 iteraciones,

**no aportan evidencia estadística suficiente** ($p = 0.080 > \alpha = 0.05$) para afirmar que la perturbación gravitacional lunar influya de forma detectable en la sismicidad global de gran magnitud.

### 3.2. Lo que el estudio sí demuestra
1. **Refuta el mito de la "Luna Llena":** el Test de Schuster sobre fase lunar no muestra agrupamiento estadísticamente significativo ($p \approx 0.45$, NB03).
2. **Cuantifica el sesgo geométrico de fondo** del sistema Tierra-Luna en aproximadamente el 60 % bajo el Criterio de Coulomb desesgado.
3. **Caracteriza la sensibilidad del resultado al declustering**, mostrando que cualquier conclusión positiva sobre triggering requiere replicar el análisis con varios procedimientos de eliminación de réplicas.

### 3.3. Lo que el estudio NO demuestra
1. **No demuestra la ausencia** del efecto. Un $p$-value alto **no es** evidencia de la nula; es ausencia de evidencia. Un catálogo más extenso, o restricción a regímenes tectónicos específicos (e.g., dorsales mediooceánicas, fallas transformantes someras), podría revelar señales locales no detectables a escala global.
2. **No invalida los hallazgos previos** en la literatura sobre subconjuntos específicos (e.g., Tanaka 2012 sobre Tohoku 2011, Cochran et al. 2004 sobre fallas reverse someras), que operan bajo condiciones físicas distintas a las del catálogo global $M \geq 7$.

### 3.4. Implicaciones para investigación futura
Cualquier publicación derivada de este trabajo debe:
- Reportar el p-value bajo **al menos dos procedimientos de declustering** independientes.
- Calibrar la línea base nula **empíricamente** (no asumir 50 %).
- Distinguir explícitamente entre la magnitud de la señal ($\Delta$ porcentaje sobre la nula) y su significancia estadística ($p$-value).
- Considerar el efecto como, en el mejor caso, un **modulador de segundo orden** sobre fallas críticamente cargadas, no como un mecanismo causal primario.

---

## 4. Reproducibilidad

Toda la cadena de cómputo está implementada en el paquete instalable [`src/lunar_trigger/`](src/lunar_trigger/) y verificada por 23 tests automatizados ([tests/](tests/)). Los resultados de la matriz de sensibilidad están en [results/statistics/declustering_sensitivity.csv](results/statistics/declustering_sensitivity.csv) y la figura comparativa en [results/figures/declustering_null_distributions.png](results/figures/declustering_null_distributions.png).

Para reproducir desde cero:

```bash
pip install -e .
jupyter nbconvert --execute --inplace notebooks/11_declustering_sensitivity.ipynb
pytest tests/
```
