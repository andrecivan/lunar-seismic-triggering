# Revisiting Lunar Tidal Triggering of Mega-Earthquakes: The Role of Aftershock Clustering and Geometric Bias in a 30-Year Global Catalog

**Iván Andrés Mena Contreras**
*ivan.andrem@gmail.com*

**Fecha del manuscrito:** 2026-04-21
**Estado:** Draft inicial — pendiente de revisión técnica y ajuste de referencias.

---

## Abstract

Re-examinamos la hipótesis del *Lunar Tidal Triggering* sobre el catálogo global de mega-terremotos $M \geq 7.0$ del USGS (1995-2024, $N = 333$). El análisis combina astrometría topocéntrica (JPL DE421 / Skyfield), proyección del tensor de marea sobre planos nodales reales bajo el Criterio de Falla de Coulomb ($\Delta\text{CFS} = \Delta\tau + \mu\,\Delta\sigma_n$, $\mu = 0.4$), y validación estadística mediante Monte Carlo *time-shuffling* con $10^3$ permutaciones. Una versión preliminar del análisis sin declustering aportaba un p-value de $0.025$, aparente rechazo de la hipótesis nula. Tres pruebas de robustez anulan esa conclusión:

1. **Declustering Gardner-Knopoff parametrizado por magnitud** (Helmstetter & Sornette, 2003) eleva el p-value monotónicamente: $0.025 \to 0.050 \to 0.080$, cruzando el umbral $\alpha = 0.05$.
2. **La línea base nula no se centra en el 50 %** sino en aproximadamente el **60 %** (catálogo global) o el **65 %** (subconjunto thrust), revelando un sesgo geométrico intrínseco del sistema Tierra-Luna que no había sido cuantificado en estudios previos.
3. **Estratificación por mecanismo focal** restringida a fallas de empuje ($N = 177$ tras declustering) — el régimen físicamente más sensible al *unclamping* lunar — produce $p = 0.268$, observado $67.23\%$ vs. media nula $65.17\%$ (exceso de $+2.06$ pp, dentro del IC 95 %).

**Conclusión:** los datos disponibles no aportan evidencia estadística suficiente para afirmar que la marea lunar module la ocurrencia de mega-terremotos a escala global, ni siquiera en zonas de subducción. La aparente correlación reportada en versiones preliminares es atribuible a (a) agrupamiento espacio-temporal de réplicas y (b) calibración incorrecta de la línea base nula. Como subproducto, este trabajo entrega una **caracterización cuantitativa del sesgo geométrico de fondo** ($\approx 60$–$65 \%$) que debe estandarizarse en futuros estudios de tidal triggering.

**Palabras clave:** tidal triggering, Coulomb failure stress, Monte Carlo, declustering, Gardner-Knopoff, sesgo geométrico, sismicidad global.

---

## 1. Introducción

La hipótesis de que las mareas terrestres y lunares pueden modular la ocurrencia de terremotos tiene más de un siglo de historia (Schuster, 1897). Aunque el efecto físico esperado es de segundo orden — la amplitud de los esfuerzos tidales crustales ($\sim 10^3$–$10^4$ Pa) es órdenes de magnitud menor que las acumulaciones tectónicas seculares ($\sim 10^6$–$10^7$ Pa) — la posibilidad de que actúe como gatillo terminal sobre fallas críticamente cargadas justifica el escrutinio recurrente.

Estudios contemporáneos han reportado señales positivas en subconjuntos específicos: Cochran et al. (2004) sobre fallas reverse someras, Tanaka (2012) sobre el segmento Tohoku 2011, Ide et al. (2016) sobre eventos $M \geq 8.2$ con catálogos extendidos. Otros trabajos no han logrado replicar significancia a escala global (Vidale et al., 1998; Hartzell & Heaton, 1989). La diversidad de resultados sugiere que (i) el efecto, si existe, es altamente sensible al régimen tectónico y al subconjunto analizado, y (ii) las metodologías estadísticas empleadas — particularmente el tratamiento de la dependencia entre eventos — son críticas para la robustez del resultado.

Este trabajo aporta tres elementos al debate:

1. **Una replicación independiente** del análisis tidal-triggering sobre el catálogo global $M \geq 7$ usando el pipeline geomecánico estándar (proyección sobre planos nodales reales, Coulomb).
2. **Una matriz de sensibilidad** del p-value frente a tres niveles progresivos de declustering, que cuantifica el grado en que las réplicas inflan la significancia aparente.
3. **Una calibración explícita de la línea base nula** mediante Monte Carlo, que revela un sesgo geométrico intrínseco del sistema Tierra-Luna no caracterizado en estudios previos y que invalida la comparación contra el 50 % asumida implícitamente en muchos análisis.

---

## 2. Datos y Métodos

### 2.1 Catálogo sísmico

Se obtuvo el catálogo global de eventos $M \geq 7.0$ del *USGS ComCat* (`https://earthquake.usgs.gov/fdsnws/event/1/query`) para el período 1995-01-01 a 2024-12-31, restringido a profundidades crustales ($\leq 70$ km), totalizando $N = 333$ eventos con mecanismo focal disponible. Para cada evento se extrajeron los dos planos nodales (`strike`, `dip`, `rake`) del producto `moment-tensor` o `focal-mechanism`.

### 2.2 Cálculo de la marea topocéntrica

La posición geocéntrica de la Luna se calculó con la efeméride DE421 (Folkner et al., 2009) vía la biblioteca Skyfield. Para cada evento, se determinó la altitud y distancia aparente de la Luna desde el hipocentro (latitud, longitud, profundidad como elevación negativa) sobre el elipsoide WGS84. Las componentes crudas del esfuerzo tidal se modelaron proporcionales a $1/r^3$ (gradiente del potencial gravitacional):

$$\sigma_{\text{raw}} \propto \frac{\sin(\text{alt})}{r^3}, \qquad \tau_{\text{raw}} \propto \frac{\cos(\text{alt})}{r^3}$$

Las constantes físicas se omiten porque solo importa la variación temporal relativa para el test estadístico.

### 2.3 Proyección sobre planos nodales y Criterio de Coulomb

Para cada plano nodal con buzamiento $\delta$, las componentes proyectadas son:

$$\sigma_n = \sigma_{\text{raw}} \cos\delta - \tau_{\text{raw}} \sin\delta$$
$$\tau = \sigma_{\text{raw}} \sin\delta + \tau_{\text{raw}} \cos\delta$$

El Coulomb Failure Stress se ensambló como $\Delta\text{CFS} = \Delta\tau + \mu\,\Delta\sigma_n$ con coeficiente de fricción estática $\mu = 0.4$ (valor estándar para fallas de subducción; sensibilidad $\mu \in [0.2, 0.6]$ no altera las conclusiones cualitativas).

**Resolución de la ambigüedad nodal.** Dado que cada mecanismo focal admite dos planos matemáticamente equivalentes (NP1, NP2) y que una elección sesgada (e.g., $\max$) infla artificialmente la fracción de eventos favorables, se adoptó el promedio simétrico:

$$\Delta\text{CFS}_{\text{unbiased}} = \frac{\Delta\text{CFS}_1 + \Delta\text{CFS}_2}{2}$$

### 2.4 Modelo nulo Monte Carlo (time-shuffling)

Para evaluar la significancia del exceso observado de eventos con $\Delta\text{CFS} > 0$, se construyó un modelo nulo barajando aleatoriamente los instantes de ocurrencia ($10^3$ permutaciones) y recalculando para cada permutación la fracción esperada bajo $H_0$. El p-value empírico (one-sided) se definió como la fracción de iteraciones cuya métrica nula iguala o supera la observada.

**Optimización computacional.** Se implementó una pre-computación matricial $N \times N$ en la que $M[i, j] = \Delta\text{CFS}$ para el evento en la ubicación $i$ ocurriendo en el tiempo $j$. La matriz se calcula con una sola invocación vectorizada a Skyfield ($N^2 = 110\,889$ evaluaciones para $N = 333$), y cada iteración Monte Carlo se reduce a una indexación por permutación, ejecutada en paralelo con `joblib`. El procedimiento es $\sim 10\times$ más rápido que el bucle ingenuo.

### 2.5 Declustering

Se aplicaron tres procedimientos progresivamente más estrictos:

- **Sin declustering** (control negativo): catálogo intacto.
- **Gardner-Knopoff fijo:** ventanas constantes de 100 km y 30 días.
- **Gardner-Knopoff variable:** ventanas $R(M)$, $T(M)$ crecientes con la magnitud, con la parametrización de Helmstetter & Sornette (2003) sobre el formalismo original de Gardner & Knopoff (1974). Esta es la opción estándar en sismicidad estadística contemporánea (Marsan & Lengliné, 2008; Zaliapin & Ben-Zion, 2013).

Implementación vectorizada con `scipy.spatial.cKDTree` para la búsqueda de vecinos espaciales (complejidad $O(N \log N)$).

### 2.6 Estratificación por mecanismo (Thrust Only)

El subconjunto thrust se definió como aquellos eventos con al menos uno de los dos rakes nodales en $[45^\circ, 135^\circ]$ (convención Aki-Richards / Frohlich 1992). Sobre este subconjunto se aplicó Gardner-Knopoff variable y Monte Carlo $10^3$ con la misma metodología.

### 2.7 Reproducibilidad

Todo el código está implementado en el paquete instalable `lunar_trigger` (Python ≥ 3.10), con 23 tests automatizados (`pytest`) que cubren: (a) paridad bit-a-bit del Coulomb desesgado contra la formulación legacy, (b) propiedades de los algoritmos de declustering (idempotencia, monotonicidad), (c) reproducibilidad por semilla del Monte Carlo. El código y los notebooks ejecutables están disponibles en el repositorio del proyecto.

---

## 3. Resultados

### 3.1 Análisis preliminar sin declustering

Sobre el catálogo completo de 333 eventos, la fracción observada de sismos con $\Delta\text{CFS}_{\text{unbiased}} > 0$ es del **64.86 %**. El modelo nulo Monte Carlo arrojó una media de **60.42 %** y un p-value empírico de **0.025** (Tabla 1). Bajo el criterio convencional $\alpha = 0.05$, este resultado constituye un aparente rechazo de la hipótesis nula.

### 3.2 Sensibilidad al declustering

La aplicación de procedimientos sucesivos de eliminación de réplicas eleva el p-value monotónicamente:

**Tabla 1.** Resultados del análisis Monte Carlo bajo tres niveles de declustering.

| Procedimiento | $N$ | % $\Delta$CFS > 0 (real) | Media nula | $\sigma_{\text{nulo}}$ | **p-value** | Decisión ($\alpha = 0.05$) |
|---|---|---|---|---|---|---|
| Sin declustering | 333 | 64.86 % | 60.42 % | 2.40 % | **0.025** | Rechaza $H_0$ |
| Gardner-Knopoff fijo (100 km / 30 d) | 309 | 63.75 % | 59.80 % | 2.46 % | **0.050** | Borderline |
| Gardner-Knopoff variable (Helmstetter-Sornette) | 288 | 63.19 % | 59.42 % | 2.58 % | **0.080** | **No rechaza $H_0$** |

El cruce del umbral $\alpha = 0.05$ ocurre exactamente al aplicar la variante fija; bajo el procedimiento parametrizado por magnitud, la significancia desaparece sin ambigüedad.

### 3.3 El sesgo geométrico del 60 %

Una observación adicional, producto colateral del análisis Monte Carlo, requiere atención metodológica explícita. La línea base bajo $H_0$ no se centra en el 50 % sugerido por un argumento ingenuo de simetría, sino en aproximadamente el **60 %** en las tres configuraciones (Tabla 1). Este desplazamiento es una propiedad geométrica intrínseca del sistema Tierra-Luna, originada en la conjunción de tres asimetrías:

1. La distribución latitudinal de las zonas de subducción que producen $M \geq 7$ se concentra en el Cinturón de Fuego del Pacífico.
2. Los planos de falla activos tienen rangos acotados de buzamiento ($20^\circ$–$60^\circ$); no son uniformes en $[0^\circ, 90^\circ]$.
3. La declinación lunar oscila entre $\pm 28.5^\circ$, privilegiando configuraciones geométricas particulares desde las zonas sísmicas activas.

**Implicación crítica:** la comparación del 64.86 % observado contra el 50 % es estadísticamente incorrecta. El exceso real es de solo $\approx 4.5$ puntos porcentuales sobre la línea base geométrica, no de $\approx 14.9$ puntos sobre el 50 % ingenuo.

### 3.4 Estratificación por mecanismo: análisis Thrust-Only

La hipótesis específica de que el efecto lunar pudiera ser detectable en zonas de subducción — donde la geometría de bajo buzamiento expone área amplia al vector lunar y la componente normal modula directamente la fricción — se evaluó sobre el subconjunto thrust del catálogo. La distribución empírica de rakes confirma que el rango $[45^\circ, 135^\circ]$ aísla apropiadamente el régimen de empuje (rake medio del catálogo: $92^\circ$, consistente con thrust puro).

**Tabla 2.** Composición por mecanismo focal del catálogo global ($N = 333$).

| Mecanismo | $N$ | Fracción |
|---|---|---|
| Thrust (rake $\in [45^\circ, 135^\circ]$) | 204 | 61.3 % |
| Strike-slip | 85 | 25.5 % |
| Normal (rake $\in [-135^\circ, -45^\circ]$) | 44 | 13.2 % |

Tras aplicar Gardner-Knopoff variable al subconjunto thrust, restan $N = 177$ mainshocks independientes (eliminación del 13.2 % por agrupamiento). El resultado del Monte Carlo $10^3$ se resume en la Tabla 3.

**Tabla 3.** Resultados del análisis estratificado por mecanismo focal (Thrust Only, GK variable, MC $10^3$).

| Métrica | Valor |
|---|---|
| Eventos analizados ($N$) | 177 |
| % observado $\Delta\text{CFS} > 0$ | 67.23 % |
| Media nula MC | 65.17 % |
| Desviación nula | 3.00 % |
| IC 95 % del nulo | [59.87 %, 71.19 %] |
| Exceso vs. nulo | +2.06 pp |
| **p-value empírico** | **0.268** |

El observado cae **dentro del intervalo de confianza al 95 % del modelo nulo**, y el exceso es menor que la desviación estándar nula (efectivamente, < $1\sigma$). Adicionalmente, la línea base geométrica sube del 60 % global al 65 % en thrust-only, confirmando que la estratificación por mecanismo amplifica — no atenúa — el sesgo geométrico de fondo.

---

## 4. Discusión

### 4.1 Convergencia de evidencia hacia el desmentido

Los tres análisis convergen en la misma dirección. (a) El test de fase lunar (no presentado en detalle aquí; ver materiales suplementarios) arroja $p \approx 0.45$ contra el Test de Schuster, descartando el mito popular de la "Luna Llena". (b) El análisis Coulomb global, una vez correctamente declusterizado, no alcanza significancia estadística. (c) La estratificación restringida al mecanismo físicamente más sensible — fallas de empuje, donde el *unclamping* normal lunar opera de forma más eficiente — no revela ninguna señal: el exceso es de $< 1\sigma$ y el p-value de $0.27$ está lejos de cualquier convención de significancia.

Esta convergencia es metodológicamente más informativa que cualquier resultado individual. Si el efecto lunar tuviera una magnitud físicamente relevante, debería emerger al menos en el régimen tectónico más favorable. Su ausencia en thrust-only constituye una **falsificación efectiva** de la hipótesis a la escala y resolución del catálogo analizado.

### 4.2 El sesgo geométrico como contribución metodológica

El hallazgo paradójicamente más sólido del estudio es la caracterización cuantitativa de la línea base geométrica nula. Que un catálogo de mega-terremotos terrestres mostrara — bajo asignación temporal puramente aleatoria — una probabilidad del $\sim 60$–$65 \%$ de configuraciones favorables al Coulomb lunar **no era resultado esperado** por análisis de simetría a primer orden, pero es robusto y reproducible al variar la semilla del Monte Carlo, los parámetros de declustering y el subconjunto analizado.

Esto tiene implicaciones directas para la literatura previa:

1. **Estudios que comparan fracciones observadas contra el 50 % deben ser revisitados.** Una fracción observada del $\sim 60 \%$ no es necesariamente evidencia de modulación tidal; puede ser exactamente lo esperado bajo $H_0$ una vez calibrada la línea base.
2. **La línea base es subset-dependiente:** estratificar por mecanismo, profundidad o región altera el sesgo geométrico. Cada análisis requiere su propia calibración Monte Carlo.
3. **Los tests paramétricos clásicos (Schuster, binomial contra 0.5) están sub-especificados** para este problema, y los reportes basados únicamente en ellos deben interpretarse con cautela.

### 4.3 Comparación con la literatura

El presente resultado negativo es consistente con los reportes de Vidale et al. (1998) y Hartzell & Heaton (1989) a escala global, y no contradice necesariamente los hallazgos positivos de Cochran et al. (2004), Tanaka (2012) e Ide et al. (2016), que operan sobre subconjuntos físicos distintos (fallas reverse someras de área limitada, segmentos pre-mainshock específicos, eventos $M \geq 8.2$). La hipótesis residual — que el efecto sea detectable solo en regímenes muy específicos no agregables a escala global — permanece abierta.

### 4.4 Limitaciones

1. **Tamaño del catálogo.** $N = 333$ ($177$ thrusts independientes) es modesto para detectar efectos de magnitud $\lesssim 5 \%$ con potencia razonable. La replicación con catálogos extendidos en magnitud ($M \geq 6.5$, $N \sim 10^3$) o en período (incluyendo ISC pre-1995) es deseable.
2. **Valores anómalos en el catálogo USGS.** Se detectó al menos un evento con `rake1 = 359°` (probable artefacto de wrap-around en la fuente). El filtrado por disyunción lógica entre los dos planos nodales mitiga el problema, pero recomendamos auditoría sistemática del catálogo de mecanismos focales en futuras versiones del estudio.
3. **Ausencia de modelado de mareas oceánicas.** El presente análisis modela exclusivamente la marea sólida directa (gravitacional). En zonas costeras, la carga oceánica puede dominar la perturbación de Coulomb, y su omisión podría enmascarar señales locales (ver Wilcock, 2001 para tratamiento detallado).
4. **Coeficiente de fricción único.** Aunque se verificó robustez para $\mu \in [0.2, 0.6]$, valores extremos o variabilidad espacial de $\mu$ podrían alterar conclusiones marginales.

### 4.5 Líneas futuras

1. Replicar la matriz de sensibilidad sobre catálogos regionales específicos (Japón, Chile, Aleutianas) donde la densidad de mecanismos focales y la calidad del catálogo son superiores.
2. Incorporar declustering no paramétrico (Zaliapin & Ben-Zion, 2013; nearest-neighbor) como tercer nivel de robustez.
3. Acoplar el modelo de marea sólida con cargas oceánicas realistas (FES2014, GOT4.10) para subconjuntos costeros.
4. Evaluar métodos bayesianos (e.g., factor de Bayes contra Schuster) que cuantifiquen evidencia *a favor* de la nula, no solo ausencia de evidencia contra ella.

---

## 5. Conclusiones

1. **No se rechaza la hipótesis nula** de ausencia de modulación lunar sobre la sismicidad global $M \geq 7$. El p-value bajo Gardner-Knopoff variable es $0.080$ (catálogo completo) y $0.268$ (subconjunto thrust).
2. La aparente significancia del análisis sin declustering ($p = 0.025$) es **artefacto del agrupamiento de réplicas**, eliminada por procedimientos estándar de declustering.
3. La línea base geométrica del modelo nulo está cerca del **60 % en el catálogo global** y del **65 % en el subconjunto thrust**, no del 50 % asumido implícitamente en muchos estudios. Comparar fracciones observadas contra el 50 % es estadísticamente incorrecto.
4. La estratificación por mecanismo de empuje — el régimen físicamente más favorable al *unclamping* lunar — no revela ninguna señal residual. Esto constituye una falsificación efectiva del efecto a la escala del catálogo.
5. El estudio aporta como subproducto una caracterización cuantitativa del sesgo geométrico de fondo y un pipeline reproducible (paquete `lunar_trigger`, 23 tests automatizados) aplicable a otros catálogos.

---

## Referencias

- Aki, K. & Richards, P. G. (2002). *Quantitative Seismology* (2nd ed.). University Science Books.
- Cochran, E. S., Vidale, J. E. & Tanaka, S. (2004). Earth tides can trigger shallow thrust fault earthquakes. *Science*, 306(5699), 1164-1166.
- Folkner, W. M., Williams, J. G. & Boggs, D. H. (2009). The Planetary and Lunar Ephemeris DE 421. *JPL Interplanetary Network Progress Report*, 42-178.
- Frohlich, C. (1992). Triangle diagrams: ternary graphs to display similarity and diversity of earthquake focal mechanisms. *Physics of the Earth and Planetary Interiors*, 75(1-3), 193-198.
- Gardner, J. K. & Knopoff, L. (1974). Is the sequence of earthquakes in southern California, with aftershocks removed, Poissonian? *Bulletin of the Seismological Society of America*, 64(5), 1363-1367.
- Hartzell, S. & Heaton, T. (1989). The fortnightly tide and the tidal triggering of earthquakes. *Bulletin of the Seismological Society of America*, 79(4), 1282-1286.
- Helmstetter, A. & Sornette, D. (2003). Foreshocks explained by cascades of triggered seismicity. *Journal of Geophysical Research*, 108(B10), 2457.
- Ide, S., Yabe, S. & Tanaka, Y. (2016). Earthquake potential revealed by tidal influence on earthquake size-frequency statistics. *Nature Geoscience*, 9(11), 834-837.
- Marsan, D. & Lengliné, O. (2008). Extending earthquakes' reach through cascading. *Science*, 319(5866), 1076-1079.
- Schuster, A. (1897). On lunar and solar periodicities of earthquakes. *Proceedings of the Royal Society of London*, 61, 455-465.
- Tanaka, S. (2012). Tidal triggering of earthquakes prior to the 2011 Tohoku-Oki earthquake ($M_w$ 9.1). *Geophysical Research Letters*, 39(7), L00G26.
- Vidale, J. E., Agnew, D. C., Johnston, M. J. S. & Oppenheimer, D. H. (1998). Absence of earthquake correlation with Earth tides: An indication of high preseismic fault stress rate. *Journal of Geophysical Research*, 103(B10), 24567-24572.
- Wilcock, W. S. D. (2001). Tidal triggering of microearthquakes on the Juan de Fuca Ridge. *Geophysical Research Letters*, 28(20), 3999-4002.
- Zaliapin, I. & Ben-Zion, Y. (2013). Earthquake clusters in southern California I: Identification and stability. *Journal of Geophysical Research: Solid Earth*, 118(6), 2847-2864.

---

## Materiales suplementarios

- Repositorio del código: paquete `lunar_trigger` con módulos `data/`, `physics/`, `stats/`, `utils/`.
- Notebooks ejecutables: `notebooks/11_declustering_sensitivity.ipynb` (matriz de declustering global), `notebooks/12_thrust_only_validation.ipynb` (estratificación por mecanismo).
- Tabla completa de KPIs: `results/statistics/declustering_sensitivity.csv`, `results/statistics/thrust_only_kpis.csv`.
- Figuras: `results/figures/declustering_null_distributions.png`, `results/figures/11_thrust_only_validation.png`.
- Tests automatizados: `tests/` (23 tests cubriendo Coulomb, declustering, Monte Carlo).
