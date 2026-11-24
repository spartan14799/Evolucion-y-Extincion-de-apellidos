# `src/core/graphTopologies.jl` — Generadores de topologías del grafo

## Qué hace

Construye el diccionario de adyacencia `compatibles::Dict{Int,Vector{Int}}` que
consume el modelo por parejas (`src/coupleModel.jl`): `compatibles[i]` es la
lista de apellidos con los que el apellido `i` puede emparejarse. Cada función
devuelve una topología distinta (completo, cadena, anillo, estrella, grid,
bipartito, tripartito, componentes).

## Cómo lo hace

Cada apellido es un vértice identificado por un entero `1..n`. Por convención,
cada apellido es compatible consigo mismo (auto-lazo) salvo que se pase
`autolazo=false`. Las topologías son simétricas en la compatibilidad (si `i`
puede con `j`, `j` puede con `i`), aunque el motor calcula las probabilidades
`P(i->j)` por fila por separado.

Nota de dinámica: con topología fija y un nodo por apellido, la deriva lleva a
fijar un único apellido si el grafo es conexo. El único equilibrio real con `k`
apellidos es un grafo con `k` componentes aisladas (`topo_componentes`); el
resto producen coexistencia transitoria (lo interesante es el orden y el tiempo
de extinción).

## Qué expone

| Función | Firma | Topología |
|:--|:--|:--|
| `topo_completo` | `topo_completo(n; autolazo=true)` | Todos compatibles con todos (baseline). |
| `topo_cadena` | `topo_cadena(n; autolazo=true)` | Camino `1-2-...-n`; modela apellidos puente. |
| `topo_anillo` | `topo_anillo(n; autolazo=true)` | Cadena cerrada (ciclo). |
| `topo_estrella` | `topo_estrella(n; centro=1, autolazo=true)` | Un centro conectado a todas las hojas. |
| `topo_grid` | `topo_grid(filas, cols; autolazo=true)` | Malla con vecinos ortogonales; id `(r-1)*cols+c`. |
| `topo_bipartito` | `topo_bipartito(grupoA, grupoB; cruzado=true, autolazo=true)` | Bipartito cruzado, o dos componentes aisladas (`cruzado=false`, caso Romeo y Julieta). |
| `topo_tripartito` | `topo_tripartito(gA, gB, gC; ciclico=true, autolazo=true)` | Tres grupos en anillo (`ciclico=true`) o en cadena A-B-C. |
| `topo_componentes` | `topo_componentes(bloques::Vector{Int}...; autolazo=true)` | Bloques completamente aislados (un superviviente por bloque). |

Función privada: `_agregar_autolazo!(compatibles)`.

## Dependencias

Ninguna. La salida alimenta `src/coupleModel.jl`.

## Uso

```julia
include(joinpath("src", "core", "graphTopologies.jl"))

topo_cadena(4)                 # Dict(1=>[1,2], 2=>[1,2,3], 3=>[2,3,4], 4=>[3,4])
topo_bipartito([1,2],[3,4]; cruzado=false)   # Romeo y Julieta: 2 componentes
```
