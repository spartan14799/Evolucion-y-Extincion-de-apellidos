# `src/experiments/daneBgwRunners.jl` — Orquestación Monte Carlo BGW (DANE)

## Qué hace

Orquesta simulaciones Monte Carlo del modelo BGW aplicado a los apellidos del
Censo DANE 2018. Separa la lógica de simulación de la presentación: alimenta el
notebook `notebooks/results/bgw_colombia_dane.jl`, que sólo dibuja. No
reimplementa el modelo; llama a `simular_generaciones` de `src/bgwModel.jl`.

## Cómo lo hace

`correr_montecarlo_colombia` repite la simulación de un linaje desde `z0`
portadores y recoge la generación de extinción de cada corrida, la probabilidad
de extinguirse en `<= gen_eval` generaciones y estadísticos del tiempo de
extinción. `calcular_medianas_10` hace lo propio para los 10 apellidos a la vez,
usando `CONTEOS_DANE[i]` como población inicial de cada uno.

## Qué expone

| Función | Firma | Devuelve |
|:--|:--|:--|
| `correr_montecarlo_colombia` | `correr_montecarlo_colombia(modelo, z0, n_sims, gen_eval, semilla)` | `NamedTuple` con `tiempos_ext`, `trayectorias_muestra` (hasta 6), `extintos_en_X`, `prob_ext_X`, `mediana_ext`, `media_ext`, `max_gen_sim`. |
| `calcular_medianas_10` | `calcular_medianas_10(modelo; n_reps=30, semilla=123)` | `Vector{Float64}` con la generación mediana de extinción de cada apellido, en el orden de `APELLIDOS_DANE`. |

## Dependencias

- Externas: `Random`, `StatsBase`.
- Del proyecto (cargados antes): `src/core/numericalMethods.jl`,
  `src/bgwModel.jl` y `src/experiments/daneData.jl`.

## Uso

```julia
include(joinpath("src", "core", "numericalMethods.jl"))
include(joinpath("src", "bgwModel.jl"))
include(joinpath("src", "experiments", "daneData.jl"))
include(joinpath("src", "experiments", "daneBgwRunners.jl"))

m = ModeloPBGW(0.5; max_hijos = 15)
res = correr_montecarlo_colombia(m, CONTEOS_DANE[1], 100, 20, 42)
res.prob_ext_X, res.mediana_ext
medianas = calcular_medianas_10(m)
```
