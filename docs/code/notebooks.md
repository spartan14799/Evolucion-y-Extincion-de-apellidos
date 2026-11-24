# `notebooks/` — Notebooks Pluto

Notebooks interactivos de Pluto. Se dividen en dos carpetas según su propósito:
`basic/` ilustra cómo funcionan los modelos; `results/` aplica los modelos a
datos reales del Censo DANE 2018.

Todos cargan el código exclusivamente de `src/` mediante `include()`: **no
redefinen los modelos**. Lo único que vive en los notebooks es la presentación
(sliders, figuras, tablas). Se abren con el lanzador `scripts/pluto.jl`.

---

## `notebooks/basic/bgw_model.jl` — Modelo BGW (ilustración)

Laboratorio interactivo del modelo de Bienaymé-Galton-Watson. Secciones: los
tres regímenes según la media `m`, la función generadora y el punto fijo `q`,
trayectorias de linajes, validación Monte Carlo, la transición de fase en
`m = 1`, y una caracterización explícita de la cadena de Markov (matriz de
transición `P`, vector de estado, flujo de probabilidad y diagrama de estados).

Carga `src/core/numericalMethods.jl` y `src/bgwModel.jl`. Usa `ModeloPBGW`,
`esperanza`, `es_extincion_segura`, `fgp`, `probabilidad_extincion`,
`simular_generaciones`, `frecuencia_extincion`. El bloque de cadena de Markov
tiene su propio slider de `λ` (`lambda_markov_slider`).

Abrir con `make bgw`.

---

## `notebooks/basic/couple_model.jl` — Modelo por parejas (ilustración)

Laboratorio interactivo del modelo por parejas. Elige una topología y sus
parámetros; construye el grafo, lo dibuja, corre una simulación y muestra su
evolución: conteos y diversidad por generación, grafo final, experimentos
agregados (muchas corridas), área apilada de proporciones, matriz de
probabilidad de emparejamiento `P(i->j)` y distribución del tiempo hasta
absorción. Incluye una guía de experimentos sugeridos (Romeo y Julieta, apellido
puente, etc.).

Carga `src/coupleModel.jl`, `src/core/graphTopologies.jl`,
`src/core/experimentMetrics.jl`, `src/core/visualization.jl`.

Abrir con `make parejas`.

---

## `notebooks/results/bgw_colombia_dane.jl` — Resultados BGW con datos DANE

Aplica el modelo BGW a los 10 apellidos más frecuentes de Colombia (Censo DANE
2018): probabilidad y tiempo de extinción por apellido, histograma de tiempos,
trayectorias en escala logarítmica y comparativa de los 10 apellidos. La
orquestación Monte Carlo vive en `src/experiments/daneBgwRunners.jl`; el notebook
sólo elige parámetros y dibuja.

Carga `src/core/numericalMethods.jl`, `src/bgwModel.jl`,
`src/experiments/daneData.jl`, `src/experiments/daneBgwRunners.jl`.

Abrir con `make bgw-dane`.

---

## `notebooks/results/couple_genealogy_dane.jl` — Árbol genealógico (datos DANE)

Dibuja el árbol genealógico del modelo por parejas sobre los 10 apellidos del
Censo DANE 2018: cada fila es una generación, cada círculo una persona (color =
apellido), y las líneas conectan padres con hijos. Muestra el apellido dominante
y la cronología de extinción. La simulación vive en
`src/experiments/daneCoupleRunners.jl` (`simular_arbol_dane`); el notebook sólo
dibuja el árbol a partir del `historial`.

Carga `src/coupleModel.jl`, `src/core/graphTopologies.jl`,
`src/experiments/daneData.jl`, `src/experiments/daneCoupleRunners.jl`.

Abrir con `make parejas-dane`.
