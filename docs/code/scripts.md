# `scripts/` — Scripts ejecutables

Tres scripts que se corren desde la raíz del proyecto con `julia --project=. ...`
(o vía `make`). Cada uno carga de `src/` lo que necesita mediante `include()`.

---

## `scripts/pluto.jl` — Lanzador de notebooks Pluto

**Qué hace.** Lanzador multiplataforma de Pluto.jl. Sin argumentos abre el
servidor Pluto (eliges el notebook en el navegador); con un argumento abre un
notebook concreto.

**Cómo lo hace.** Usa el intérprete con el que se lanza (no fija ruta de Julia),
ejecuta `Pkg.instantiate()` (idempotente) y llama a `Pluto.run()`. Mapea claves a
rutas en el `Dict` `NOTEBOOKS`.

**Claves disponibles:**

| Clave | Notebook | Make |
|:--|:--|:--|
| `bgw` | `notebooks/basic/bgw_model.jl` | `make bgw` |
| `parejas` | `notebooks/basic/couple_model.jl` | `make parejas` |
| `bgw-dane` | `notebooks/results/bgw_colombia_dane.jl` | `make bgw-dane` |
| `parejas-dane` | `notebooks/results/couple_genealogy_dane.jl` | `make parejas-dane` |

**Dependencias.** `Pkg`, `Pluto`. No incluye `src/`.

**Uso.** `julia --project=. scripts/pluto.jl bgw`

---

## `scripts/run_simulation.jl` — Demo BGW por consola

**Qué hace.** Simulación de demostración del modelo BGW: para un modelo de
natalidad alta (`Poisson λ=1.5`, supercrítico) y uno baja (`λ=0.8`,
subcrítico) imprime la media `m`, el régimen, la `q` teórica, unas trayectorias
de ejemplo y una validación Monte Carlo frente a `q`.

**Cómo lo hace.** `resumen_modelo(nombre, modelo)` imprime `esperanza` y
`probabilidad_extincion`; `main()` arma los dos modelos, genera trayectorias con
`simular_generaciones` y contrasta con `frecuencia_extincion`. Llama a `main()`
al final.

**Dependencias.** Externas: `Random`, `Printf`. Incluye
`src/core/numericalMethods.jl` y `src/bgwModel.jl`.

**Uso.** `julia --project=. scripts/run_simulation.jl` (o `make simulation`).

---

## `scripts/experimentos_parejas.jl` — Catálogo reproducible de experimentos

**Qué hace.** Corre sin Pluto un catálogo fijo de experimentos del modelo por
parejas con semillas fijas y escribe a `output/`: por cada experimento, un PNG
del grafo inicial, del conteo por apellido, de la diversidad y de la distribución
de supervivientes; más un `resumen_experimentos.csv` con el agregado de las
repeticiones.

**Cómo lo hace.** `construir_catalogo(parejas)` define 10 experimentos (completo,
romeo_julieta, cadena_puente, estrella, anillo, grid, bipartito_cruzado,
tripartito_ciclico, tripartito_cadena, componentes). `main()` corre cada uno con
`correr_repeticiones`, genera las figuras con `src/core/visualization.jl` y
escribe el CSV con `DelimitedFiles`. Parametrizable por variables de entorno.

**Variables de entorno** (defaults reproducibles):

| Variable | Controla | Default |
|:--|:--|:--|
| `EXP_N_REPS` | repeticiones por experimento | 200 |
| `EXP_PAREJAS` | parejas por apellido en el estado inicial | 3 |
| `EXP_SEMILLA_BASE` | semilla base de las repeticiones | 1 |
| `EXP_MAX_GEN` | tope de generaciones por corrida | 1000 |
| `EXP_SEMILLA_FIG` | semilla de la corrida de ejemplo (figuras) | 7 |
| `EXP_OUTDIR` | carpeta de salida | `output` |

**Dependencias.** Externas: `Random`, `StatsBase`, `Printf`, `DelimitedFiles`
(stdlib). Incluye `src/coupleModel.jl`, `src/core/graphTopologies.jl`,
`src/core/experimentMetrics.jl`, `src/core/visualization.jl`.

**Uso.** `julia --project=. scripts/experimentos_parejas.jl` (o `make experimentos`).
Ejemplo con otra configuración: `EXP_N_REPS=500 EXP_PAREJAS=5 julia --project=. scripts/experimentos_parejas.jl`.
