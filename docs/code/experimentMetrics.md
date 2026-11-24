# `src/core/experimentMetrics.jl` — Métricas y runner del modelo por parejas

## Qué hace

Post-procesa los resultados del modelo por parejas y orquesta corridas
repetidas. Dos grupos: **métricas** que resumen una corrida (diversidad, orden
de extinción, mutaciones) y un **runner** que ejecuta muchas corridas con
semillas reproducibles y agrega estadísticas. No toca el motor del modelo.

## Cómo lo hace

Todo opera sobre el `historial` y `forzadas_por_gen` de un `ResultadoParejas`.
El runner ejecuta `correr_experimento` con semillas `semilla_base .+ (0:n_reps-1)`
para reproducibilidad, resume cada corrida y agrega: probabilidad de absorción,
distribución del número de supervivientes, probabilidad de ganar por apellido, y
(opcional) con qué frecuencia un apellido "puente" se extingue primero.

## Qué expone

### Métricas

| Función | Firma | Propósito |
|:--|:--|:--|
| `apellidos_vivos_por_gen` | `apellidos_vivos_por_gen(historial)` | Nº de apellidos distintos vivos por generación (serie no creciente). |
| `conteos_por_gen` | `conteos_por_gen(historial, apellidos) -> Matrix{Int}` | Conteo de cada apellido por generación; `apellidos` fija las columnas. |
| `generacion_de_extincion` | `generacion_de_extincion(historial) -> Dict{Int,Int}` | Generación (base 0) en que cada apellido inicial llega a cero (`typemax(Int)` si no). |
| `orden_de_extincion` | `orden_de_extincion(historial)` | Apellidos ordenados por generación de extinción. |
| `total_mutaciones` | `total_mutaciones(resultado::ResultadoParejas) -> Int` | Total de emparejamientos forzados. |
| `mutaciones_por_gen` | `mutaciones_por_gen(resultado) -> Vector{Int}` | Emparejamientos forzados por generación. |
| `resumen_experimento` | `resumen_experimento(resultado; puente=nothing) -> NamedTuple` | Resumen de una corrida. |

### Runner

| Función | Firma | Propósito |
|:--|:--|:--|
| `poblacion_desde_parejas` | `poblacion_desde_parejas(parejas_por_apellido::Dict{Int,Int}) -> Vector{Int}` | Estado inicial desde el nº de parejas por apellido (población par). |
| `correr_experimento` | `correr_experimento(estado_inicial, compatibles; modo_prob=:uniforme_aristas, max_generaciones=1000, semilla=nothing)` | Una sola corrida (atajo sobre `simular_parejas`). |
| `correr_repeticiones` | `correr_repeticiones(estado_inicial, compatibles; modo_prob=:uniforme_aristas, max_generaciones=1000, n_reps=100, semilla_base=1, puente=nothing) -> NamedTuple` | N corridas con estadísticas agregadas. |

`correr_repeticiones` devuelve, entre otros: `prob_absorcion`, `prob_transitorio`,
`gen_media`, `gen_std`, `dist_n_vivos`, `prob_ganar`, `mutaciones_media`,
`prob_puente_primero`, `resumenes`.

## Dependencias

- Externas: `Random`, `StatsBase`.
- Del proyecto: `src/coupleModel.jl` (usa `ResultadoParejas`, `ParametrosParejas`,
  `simular_parejas`), que debe estar cargado antes.

## Uso

```julia
include(joinpath("src", "coupleModel.jl"))
include(joinpath("src", "core", "graphTopologies.jl"))
include(joinpath("src", "core", "experimentMetrics.jl"))

estado = poblacion_desde_parejas(Dict(1 => 3, 2 => 3))   # 2 apellidos, 3 parejas c/u
compat = topo_completo(2)
rep = correr_repeticiones(estado, compat; n_reps = 200)
rep.prob_absorcion, rep.dist_n_vivos
```
