# Revisiting Lunar Tidal Triggering: The Role of Clustering and Geometric Bias
**Autor:** Ing. Iván Andrés Mena Contreras

## 1. Resumen Ejecutivo (Abstract)
Este estudio re-evalúa la hipótesis del desencadenamiento sísmico inducido por la gravedad lunar (*Lunar Tidal Triggering*) sobre un catálogo global de mega-terremotos ($M \geq 7.0$, 1995-2024). Combina astrometría de alta precisión (JPL Horizons / Skyfield), modelado vectorial del tensor de esfuerzos sobre planos nodales reales (Criterio de Falla de Coulomb) y validación estadística mediante Monte Carlo de tipo *time-shuffling*. Una versión preliminar del análisis sin declustering aportaba un **p-value de 0.025** (rechazo aparente de la hipótesis nula). Sin embargo, dos pruebas de robustez clave anulan esa conclusión:

1. **Efecto del declustering:** al eliminar réplicas mediante el método de Gardner-Knopoff parametrizado por magnitud (Helmstetter & Sornette 2003), el p-value **crece monotónicamente hasta 0.080**, cruzando el umbral $\alpha = 0.05$ y eliminando la significancia estadística.
2. **Sesgo geométrico de fondo:** la línea base bajo la hipótesis nula no se centra en el 50 % esperado por simetría, sino en aproximadamente el **60 %**, debido a la asimetría conjunta de la distribución de fallas terrestres y la órbita lunar. La aparente "señal" del 64.86 % es solo $\approx 4.5$ puntos sobre esta línea geométrica.

**Conclusión:** *no se puede rechazar la hipótesis nula*. La aparente correlación lunar-sísmica observada en versiones preliminares es **artefacto de agrupamiento de réplicas y sesgo geométrico no calibrado**, no evidencia de un mecanismo gatillante real. El veredicto detallado se documenta en [`CONCLUSIONS_SUMMARY.md`](CONCLUSIONS_SUMMARY.md).

## 2. Arquitectura de Datos y Pipeline
El proyecto estructura un flujo de procesamiento automatizado, riguroso y reproducible, encapsulado en el paquete instalable [`src/lunar_trigger/`](src/lunar_trigger/) y validado por 23 tests automatizados ([`tests/`](tests/)):

- **Extracción de Catálogos:** consultas REST a la API del *United States Geological Survey* (USGS ComCat) con sesión HTTP reutilizable y caché en disco ([`data/usgs.py`](src/lunar_trigger/data/usgs.py)).
- **Astrometría y Efemérides:** posición topocentrica de la Luna en el instante y ubicación exacta de cada hipocentro mediante **JPL Horizons (NASA)** y **Skyfield** (elipsoide WGS84, efeméride DE421). Implementación vectorizada en [`data/ephemeris.py`](src/lunar_trigger/data/ephemeris.py).
- **Mecanismos Focales:** recuperación de los planos nodales (*Strike*, *Dip*, *Rake*) desde el catálogo *Moment Tensor* del USGS.
- **Proyección Tensorial:** cálculo de las componentes normal y cortante de la marea sobre los planos de falla, y ensamblado del Coulomb Failure Stress en [`physics/coulomb.py`](src/lunar_trigger/physics/coulomb.py).
- **Estadística:** Monte Carlo time-shuffling con pre-cómputo matricial $N \times N$ y paralelización joblib ([`stats/monte_carlo.py`](src/lunar_trigger/stats/monte_carlo.py)); declustering Gardner-Knopoff fijo y variable ([`stats/declustering.py`](src/lunar_trigger/stats/declustering.py)).

## 3. Fase 1: Desmontando el Mito de la "Luna Llena"
Inicialmente, se sometió a prueba la creencia popular de que los sismos se agrupan en sicigias (Luna Llena y Luna Nueva). Utilizando la fase lunar y el **Test de Schuster** (Notebook 03), la distribución de los sismos a nivel global arrojó un valor-$p > 0.05$ (específicamente $p \approx 0.45$): **distribución estadísticamente indiferenciable del azar absoluto**.

**Conclusión de la Fase 1:** la fase lunar (variable escalar bidimensional) es analíticamente insuficiente para probar o refutar transferencia de estrés sobre la litósfera. Es necesario un marco vectorial tridimensional.

## 4. Fase 2: Modelado Físico y el Criterio de Coulomb
Se abandonaron las variables escalares para implementar un marco geomecánico newtoniano. Dado que la amplitud de la marea gravitacional decae con el cubo de la distancia, el análisis se reconfiguró bajo un modelado vectorial proporcional a $1/r^3$ (Notebook 06).

La falla de los materiales bajo estrés crustal fue evaluada usando el **Criterio de Falla de Coulomb**:

$$
\Delta \text{CFS} = \Delta \tau + \mu \cdot \Delta \sigma_n
$$

donde $\tau$ es la perturbación cortante, $\sigma_n$ el esfuerzo normal (descompresión / *unclamping*) y $\mu$ el coeficiente de fricción estática (valor estándar $\mu = 0.4$ para fallas de subducción; sensibilidad evaluada en NB08).

## 5. Fase 3: Proyección en Planos de Falla Reales y Resolución de la Ambigüedad Nodal
Cada mecanismo focal admite **dos planos nodales matemáticamente equivalentes** (NP1 y NP2). Sin información adicional (afterslip, geología local), no es posible elegir cuál es el plano de ruptura real. Esta ambigüedad introduce un riesgo grave de sesgo de selección:

- Tomar $\max(\Delta\text{CFS}_1, \Delta\text{CFS}_2)$ es **estadísticamente sesgado**: como la marea oscila, siempre habrá un plano "favorable" por azar, lo que infla artificialmente la fracción de eventos favorecidos.
- La opción conservadora es **promediar** los dos planos: $(\Delta\text{CFS}_1 + \Delta\text{CFS}_2)/2$ (Notebook 09).

Bajo esta resolución desesgada, sin declustering, la fracción de eventos con $\Delta\text{CFS} > 0$ resulta **64.86 %**.

## 6. Fase 4: Calibración Empírica de la Línea Base — El Sesgo Geométrico del 60 %
Para juzgar si el 64.86 % constituye exceso real, se construyó un modelo nulo mediante **Monte Carlo time-shuffling** (Notebook 09): los 333 instantes de ocurrencia se barajaron aleatoriamente, se recalcularon las mareas topocentricas y se midió la fracción esperada bajo H₀.

El resultado fue inesperado y metodológicamente importante: la línea base nula **no se centra en el 50 %** que un argumento ingenuo de simetría predeciría, sino en aproximadamente el **60 %**. Este desplazamiento es una propiedad **geométrica intrínseca** del sistema Tierra-Luna, producto de la asimetría conjunta entre:

1. La distribución latitudinal de las zonas de subducción que producen $M \geq 7$.
2. Los rangos acotados de buzamiento ($dip$) de las fallas activas.
3. La oscilación de la declinación lunar ($\pm 28.5°$).

**Implicación crítica:** comparar el 64.86 % observado contra el 50 % es **estadísticamente incorrecto** y produce conclusiones falsamente positivas. El exceso real es solo $\approx 4.5$ puntos porcentuales sobre la línea base geométrica de 60.39 %, no $\approx 14.9$ puntos sobre el 50 % ingenuo. Versiones preliminares del manuscrito reportaron un p-value de **0.025** sobre esta base.

## 7. Fase 5: Análisis de Robustez — La Señal Muere bajo Declustering
El modelo nulo de time-shuffling asume **independencia temporal** entre los eventos del catálogo. Esa asunción se viola si el catálogo contiene replicas (mainshock + aftershocks): pares de sismos ligados por una causa física común (transferencia de Coulomb estática del mainshock) producirán una falsa apariencia de modulación cuando lo que existe es agrupamiento intrínseco.

Para evaluar la robustez del p-value frente a esta violación, se aplicaron tres niveles progresivamente más estrictos de declustering (Notebook 11):

| Procedimiento | N | % $\Delta$CFS > 0 | Media nula | **p-value** |
|---|---|---|---|---|
| Sin declustering (NB09) | 333 | 64.86 % | 60.39 % | **0.025** |
| Gardner-Knopoff fijo, 100 km / 30 d (NB10) | 309 | 63.75 % | 59.65 % | **0.050** |
| Gardner-Knopoff variable, Helmstetter & Sornette 2003 (NB11) | 288 | 63.19 % | 59.42 % | **0.080** |

El p-value **crece monotónicamente** a medida que el declustering se vuelve más severo, cruzando el umbral de significancia $\alpha = 0.05$ exactamente al aplicar Gardner-Knopoff fijo y alejándose definitivamente bajo el procedimiento parametrizado por magnitud. Esto demuestra que la aparente significancia inicial era **artefacto del agrupamiento de réplicas**, no evidencia de un mecanismo gatillante.

## 8. Conclusiones Finales

> **Veredicto:** *No se puede rechazar la hipótesis nula.* La marea gravitacional lunar **no aporta evidencia estadística suficiente** ($p = 0.080 > \alpha = 0.05$) para afirmar que modula la ocurrencia de mega-terremotos $M \geq 7$ en el catálogo global 1995-2024.

El estudio refuta de forma categórica el mito popular de la "Luna Llena" como gatillo sísmico, y demuestra adicionalmente que **la aparente correlación reportada en versiones preliminares es atribuible a dos factores no relacionados con la marea**:

1. **Sesgo geométrico del 60 %** en la línea base nula (versus el 50 % ingenuamente esperado), producto de la asimetría tectónica-orbital del sistema Tierra-Luna.
2. **Agrupamiento de réplicas** que viola la independencia temporal asumida por el modelo nulo.

Lo que el estudio sí establece de forma sólida es:
- Una **caracterización cuantitativa del sesgo geométrico de fondo** ($\approx 60 \%$), de utilidad metodológica para futuros estudios de tidal triggering.
- Un **pipeline reproducible y testeado** (paquete `lunar_trigger`, 23 tests pytest, notebooks ejecutables) que cualquier investigador puede aplicar a otros catálogos o regiones específicas.
- Una **demostración explícita** de que el resultado depende críticamente del procedimiento de declustering, lo cual debe ser estándar de reporte en cualquier publicación futura sobre el tema.

El detalle metodológico y las implicaciones para investigación futura se documentan en [`CONCLUSIONS_SUMMARY.md`](CONCLUSIONS_SUMMARY.md).

## 9. Reproducibilidad

```bash
# Instalación
python -m venv venv
source venv/Scripts/activate    # Windows: venv\Scripts\activate
pip install -e .[notebook,test]

# Re-ejecutar el análisis de sensibilidad
jupyter nbconvert --execute --inplace notebooks/11_declustering_sensitivity.ipynb

# Verificar la integridad del paquete
pytest tests/ -v
```

Los resultados intermedios se almacenan en [`data/processed/`](data/processed/) y los productos publicables en [`results/`](results/).
