# `src/experiments/daneData.jl` — Datos del Censo DANE 2018

## Qué hace

Archivo de **sólo datos**: los 10 apellidos más frecuentes de Colombia según el
Censo Nacional de Población y Vivienda (DANE 2018), sus conteos y una paleta de
colores. No depende de ningún modelo. Lo consumen los notebooks de
`notebooks/results/` y los dos archivos de orquestación (`daneBgwRunners.jl` y
`daneCoupleRunners.jl`).

## Cómo lo hace

Expone los datos como constantes (`const`) en orden de frecuencia descendente,
más dos helpers derivados. Hay dos versiones de los nombres: con acentos para
tablas Markdown y sin acentos (ASCII) para figuras cuya fuente no renderiza bien
los diacríticos.

## Qué expone

| Símbolo | Tipo / Firma | Contenido |
|:--|:--|:--|
| `APELLIDOS_DANE` | `const Vector{String}` | 10 apellidos con acentos, orden de frecuencia. |
| `APELLIDOS_DANE_ASCII` | `const Vector{String}` | Los mismos sin acentos. |
| `CONTEOS_DANE` | `const Vector{Int}` | Portadores de cada apellido (mismo orden). |
| `COLORES_DANE` | `const Vector{String}` | 10 colores (uno por apellido) para figuras. |
| `pesos_dane` | `pesos_dane() -> Vector{Float64}` | Conteos normalizados a suma 1. |
| `dane_dict` | `dane_dict() -> Dict{String,Int}` | Diccionario apellido => conteo (con acentos). |

## Dependencias

- Externas: `StatsBase`.
- Del proyecto: ninguna.

## Uso

```julia
include(joinpath("src", "experiments", "daneData.jl"))

APELLIDOS_DANE[1]      # "Rodríguez"
CONTEOS_DANE[1]        # 754991
pesos_dane()           # pesos relativos que suman 1
```
