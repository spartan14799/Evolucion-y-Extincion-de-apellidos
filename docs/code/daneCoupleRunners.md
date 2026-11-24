# `src/experiments/daneCoupleRunners.jl` — Orquestación del árbol por parejas (DANE)

## Qué hace

Orquesta la simulación del árbol genealógico del modelo por parejas aplicado a
los 10 apellidos del Censo DANE 2018. Corre la simulación con el motor
`siguiente_generacion_parejas` de `src/coupleModel.jl` y devuelve el historial y
los metadatos necesarios para que el notebook
`notebooks/results/couple_genealogy_dane.jl` sólo dibuje. No reimplementa el
modelo.

## Cómo lo hace

`simular_arbol_dane` construye una población inicial sobre los 10 apellidos
(equitativa o proporcional al Censo), usa un grafo completo entre ellos
(`topo_completo(10)`) y avanza generación a generación guardando cada estado.
Registra en qué generación cae cada apellido y detecta la fijación (un único
apellido vivo).

## Qué expone

| Función | Firma | Devuelve |
|:--|:--|:--|
| `simular_arbol_dane` | `simular_arbol_dane(; modo_inicial=:equitativo, parejas=10, modo_prob=:ponderado, semilla=7, max_generaciones=25)` | `NamedTuple` con `historial`, `generaciones`, `caidos` (Dict apellido => generación de extinción), `ganador` (id dominante), `fijado` (Bool), `n_total` (= `2*parejas`). |

Modos de población inicial: `:equitativo` reparte las parejas cíclicamente entre
los 10 apellidos; `:dane` muestrea el apellido de cada pareja proporcional a
`pesos_dane()`.

## Dependencias

- Externas: `Random`, `StatsBase`.
- Del proyecto (cargados antes): `src/coupleModel.jl`,
  `src/core/graphTopologies.jl` y `src/experiments/daneData.jl`.

## Uso

```julia
include(joinpath("src", "coupleModel.jl"))
include(joinpath("src", "core", "graphTopologies.jl"))
include(joinpath("src", "experiments", "daneData.jl"))
include(joinpath("src", "experiments", "daneCoupleRunners.jl"))

arbol = simular_arbol_dane(; modo_inicial = :dane, parejas = 10, semilla = 7)
APELLIDOS_DANE_ASCII[arbol.ganador], arbol.generaciones, arbol.fijado
```
