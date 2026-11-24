# Modelo por parejas

Modelo de evolución y extinción de apellidos sobre un grafo. Esta nota describe
las reglas, los modos de probabilidad, el criterio exacto de parada, el mapeo al
código y el catálogo de experimentos con sus hipótesis.

El motor vive en [`src/coupleModel.jl`](../src/coupleModel.jl); las
topologías en [`src/core/graphTopologies.jl`](../src/core/graphTopologies.jl); las
métricas y el runner en [`src/core/experimentMetrics.jl`](../src/core/experimentMetrics.jl); la
graficación en [`src/core/visualization.jl`](../src/core/visualization.jl). El catálogo
reproducible a disco está en
[`scripts/experimentos_parejas.jl`](../scripts/experimentos_parejas.jl).

## 1. El modelo

Cada **apellido es un vértice** de un grafo dirigido. Las aristas, guardadas en
`compatibles::Dict{Int,Vector{Int}}`, indican qué apellidos pueden emparejarse:
`compatibles[i]` es la lista de apellidos con los que `i` puede formar pareja
(incluido, normalmente, él mismo por el auto-lazo).

La población se parametriza por el **número de parejas** `P`: hay `N = 2P`
personas, de modo que `N` siempre es par. El estado es un `Vector{Int}` de
longitud `N` donde cada entrada es el apellido de una persona.

### Reglas de una generación (de `2P` padres a `2P` hijos)

1. **Emparejamiento sin reemplazo.** Se barajan las personas y, en orden
   aleatorio, cada persona suelta busca pareja entre las personas **todavía
   sueltas** cuyo apellido sea compatible con el suyo. La pareja concreta se
   sortea según el modo de probabilidad (ver §2), con pesos calculados **una vez
   al inicio de la generación**. Al emparejarse, las dos personas salen de la
   bolsa de sueltos.
2. **Herencia.** Por cada pareja se elige al azar (50/50) a uno de los dos
   miembros; **los dos hijos** de esa pareja toman ese apellido. El otro apellido
   no se transmite en esa pareja.
3. **Tamaño constante.** Cada pareja produce exactamente 2 hijos, así que
   `N = 2P` se conserva generación a generación.

### Mutación definida (emparejamiento forzado)

Si al terminar el emparejamiento quedan personas sueltas **sin ningún compatible
suelto disponible**, se emparejan entre sí aunque sean incompatibles. Cada uno
de esos eventos es una **mutación definida** y se registra como un par de
apellidos `(a, b)`. Son la vía por la que un apellido puede aparecer en una
pareja donde "no debería" según las aristas, y se contabilizan por generación.

## 2. Modos de probabilidad

El campo `modo_prob` de `ParametrosParejas` elige cómo se reparte la
probabilidad de que una persona de apellido `i` se empareje con cada apellido
compatible `j`. Solo se consideran apellidos **vivos** (presentes en la
generación actual): si un apellido se extingue, su arista desaparece y la
probabilidad se reparte entre las aristas restantes.

- **`:uniforme_aristas`** (por defecto): todas las aristas vivas de `i` pesan
  igual. Si `i` tiene `d` vecinos vivos, cada uno recibe peso `1/d`. Es el modo
  "neutral por estructura": la topología manda, no los tamaños.
- **`:ponderado`**: el peso de `j` es proporcional al **conteo** de `j` entre
  los compatibles de `i` (calculado una vez por generación). Favorece
  emparejarse con los apellidos más numerosos.

La función que calcula los pesos es `pesos_emparejamiento`; devuelve los
`destinos` posibles y sus `pesos` sin normalizar.

## 3. Criterio exacto de parada

La diversidad (número de apellidos distintos) es una serie **no creciente**:
solo puede bajar cuando se forma una **pareja mixta** (dos personas de apellidos
distintos y compatibles), porque ahí un apellido le gana el relevo al otro en el
volado de la herencia. Si dos personas del mismo apellido se emparejan, sus hijos
conservan ese apellido y el grupo no pierde tamaño.

De ahí el criterio exacto, implementado en `hay_parejas_mixtas_posibles`:

> Si entre los apellidos vivos **no** existe ninguna arista `i -> j` con
> `i != j`, ya no pueden formarse parejas mixtas y la diversidad no puede bajar
> más: el estado es un **equilibrio permanente**.

`simular_parejas` se detiene al ocurrir lo primero de:

1. **Fijación**: queda un único apellido. `absorbido = true`, `transitorio = false`.
2. **Equilibrio permanente**: ya no hay parejas mixtas posibles entre los
   apellidos vivos (p. ej. Romeo y Julieta con dos grupos aislados).
   `absorbido = true`, `transitorio = false`.
3. **Corte por tope**: se alcanza `max_generaciones` sin cumplir 1 ni 2. Si a
   esa altura todavía hay parejas mixtas posibles, el estado se reporta como
   `transitorio = true` (la diversidad aún podría bajar con más generaciones).

### Fijación / equilibrio / transitorio

- **Fijación**: estado absorbente con **un** apellido. Es el destino típico de
  una topología conexa con población pequeña.
- **Equilibrio** (permanente): estado absorbente con **k ≥ 2** apellidos porque
  la estructura ya no admite parejas mixtas. El único modo de garantizarlo es
  tener `k` componentes aisladas (una por apellido superviviente).
- **Transitorio**: la corrida se cortó por el tope y el estado **aún admite
  parejas mixtas**; la coexistencia observada no es estable, solo todavía no
  colapsó. `absorbido` y `transitorio` son excluyentes.

## 4. Mapeo al código

| Concepto | Símbolo en el código |
|:--|:--|
| Parámetros de la simulación | `ParametrosParejas(max_generaciones, modo_prob)` |
| Una transición de generación | `siguiente_generacion_parejas(poblacion, compatibles, rng; modo_prob)` → `(nueva_poblacion, parejas_forzadas)` |
| Simulación completa | `simular_parejas(estado_inicial, compatibles; parametros, rng)` → `ResultadoParejas` |
| Pesos de emparejamiento | `pesos_emparejamiento(i, compatibles, conteos, vivos, modo)` → `(destinos, pesos)` |
| Criterio de parada | `hay_parejas_mixtas_posibles(apellidos_vivos, compatibles)` |
| Resultado de una corrida | `ResultadoParejas` |
| Estado inicial desde nº de parejas | `poblacion_desde_parejas(Dict(apellido => n_parejas))` |

### `ResultadoParejas`

| Campo | Significado |
|:--|:--|
| `estado_final` | población en la última generación |
| `generaciones` | nº de transiciones simuladas |
| `historial` | estado de cada generación (incluida la 0) |
| `forzadas_por_gen` | parejas forzadas (mutaciones) de cada transición |
| `absorbido` | alcanzó fijación o equilibrio permanente |
| `apellidos_final` | apellidos vivos al final |
| `transitorio` | se cortó por tope con parejas mixtas aún posibles |

### Generadores de topología (`src/core/graphTopologies.jl`)

`topo_completo(n)`, `topo_cadena(n)`, `topo_anillo(n)`, `topo_estrella(n;
centro)`, `topo_grid(filas, cols)`, `topo_bipartito(gA, gB; cruzado)`,
`topo_tripartito(gA, gB, gC; ciclico)`, `topo_componentes(bloques...)`. Todos
devuelven el `Dict{Int,Vector{Int}}` de adyacencia y por defecto incluyen el
auto-lazo `i -> i`.

### Métricas (`src/core/experimentMetrics.jl`)

| Métrica | Función |
|:--|:--|
| Apellidos vivos por generación | `apellidos_vivos_por_gen(historial)` |
| Conteo de cada apellido por generación | `conteos_por_gen(historial, apellidos)` |
| Generación de extinción de cada apellido | `generacion_de_extincion(historial)` |
| Orden de extinción | `orden_de_extincion(historial)` |
| Total de mutaciones / por generación | `total_mutaciones(res)` / `mutaciones_por_gen(res)` |
| Resumen de una corrida | `resumen_experimento(res; puente)` |
| Agregado de muchas corridas | `correr_repeticiones(estado, compatibles; ...)` |

`correr_repeticiones` usa semillas reproducibles
(`semilla_base .+ (0:n_reps-1)`) y agrega: `prob_absorcion`, `prob_transitorio`,
`gen_media`/`gen_std`, `dist_n_vivos`, `prob_ganar` (por apellido),
`mutaciones_media` y `prob_puente_primero`.

## 5. Catálogo de experimentos e hipótesis

El script [`scripts/experimentos_parejas.jl`](../scripts/experimentos_parejas.jl)
corre estos experimentos con semillas fijas y guarda figuras + CSV en `output/`.

| Experimento | Topología | Hipótesis |
|:--|:--|:--|
| **completo** | `topo_completo(4)` | Baseline conexo: fija **un** apellido; cada uno gana con probabilidad ≈ su frecuencia inicial. |
| **Romeo y Julieta** | `topo_bipartito([1,2],[3,4]; cruzado=false)` | Dos componentes aisladas: equilibrio permanente de **2** apellidos de entrada, sin mutaciones. |
| **cadena / puente** | `topo_cadena(4)`, puente = 2 | Los apellidos 2 y 3 son puentes; a menudo un puente se extingue y quedan dos extremos no adyacentes coexistiendo. Aparecen mutaciones cuando un puente cae. |
| **estrella** | `topo_estrella(5; centro=1)`, puente = 1 | El centro está sobre-expuesto a parejas mixtas; cuando sobrevive tiende a dominar, pero es candidato a extinguirse primero. Hojas incompatibles entre sí fuerzan muchas mutaciones. |
| **anillo** | `topo_anillo(4)` | Como la cadena pero sin extremos; sirve para aislar el efecto de los bordes sobre el tiempo y el orden de extinción. |
| **grid** | `topo_grid(2,2)` | Vecindad local ortogonal; compara el tiempo de fijación contra el grafo completo del mismo nº de apellidos. |
| **bipartito cruzado** | `topo_bipartito([1,2],[3,4]; cruzado=true)` | Cada grupo solo se empareja con el otro (una sola componente conexa); fija un apellido, con mutaciones cuando un lado se queda sin compatibles. |
| **tripartito cíclico** | `topo_tripartito([1],[2],[3]; ciclico=true)` | Tres apellidos en anillo de grupos (una componente); fija uno. |
| **tripartito en cadena** | `topo_tripartito([1],[2],[3]; ciclico=false)`, puente = 2 | El apellido 2 es el único puente entre 1 y 3; su extinción separa el grafo. |
| **componentes** | `topo_componentes([1,2],[3,4],[5,6])` | Tres componentes aisladas: equilibrio permanente de **3** linajes, sin mutaciones. |

El CSV `output/resumen_experimentos.csv` tiene una fila por experimento con las
columnas: `experimento`, `descripcion`, `n_apellidos`, `n_reps`,
`parejas_por_apellido`, `poblacion`, `modo_prob`, `max_gen`, `semilla_base`,
`prob_absorcion`, `prob_transitorio`, `gen_media`, `gen_std`, `n_vivos_medio`,
`mutaciones_media`, `prob_puente_primero`, `dist_n_vivos`, `prob_ganar`.

## 6. Por qué este modelo conserva los grupos

La propiedad central del modelo por parejas es que la **pareja** es la unidad
reproductiva, no el individuo ni un cupo global. Como cada pareja aporta
exactamente dos hijos que heredan el apellido de uno de sus dos miembros, un
grupo cerrado (un conjunto de apellidos sin aristas hacia fuera) **conserva su
tamaño** generación a generación: las personas de ese grupo solo pueden
emparejarse entre sí, y sus hijos siguen dentro del grupo.

De esa conservación se sigue el resultado que distingue a este modelo: el
equilibrio con `k > 1` apellidos es posible y **estable** cuando el grafo tiene
`k` componentes aisladas (una por linaje superviviente). Dentro de cada
componente conexa la diversidad sí colapsa hasta fijar un apellido, pero entre
componentes nunca se pierde diversidad porque no hay parejas mixtas posibles
(ver §3). Por eso topologías como *Romeo y Julieta* (§5) mantienen dos apellidos
para siempre en lugar de derivar hacia uno solo.
