# `src/core/visualization.jl` — Utilidades de graficación

## Qué hace

Convierte las salidas del modelo por parejas en figuras: dibuja el grafo de
apellidos y grafica conteos, diversidad y estadísticas de supervivencia. No
contiene lógica de modelo; sólo presentación.

## Cómo lo hace

Usa `Plots.jl` para las gráficas y `Graphs.jl` + `GraphRecipes.jl` para dibujar
el grafo dirigido de apellidos. Las funciones de barras (`graficar_dist_n_vivos`,
`graficar_prob_ganar`) consumen directamente el `NamedTuple` que devuelve
`correr_repeticiones`.

## Qué expone

| Función | Firma | Propósito |
|:--|:--|:--|
| `dibujar_grafo` | `dibujar_grafo(compatibles; conteos=nothing, titulo="Grafo de apellidos", mostrar_autolazo=false)` | Grafo dirigido; tamaño/etiqueta de nodo según `conteos`. Maneja el caso de 0-1 nodos o sin aristas con un marcador en vez de fallar. |
| `graficar_conteos` | `graficar_conteos(historial; apellidos=nothing, titulo="Conteo por apellido")` | Una curva por apellido con su conteo por generación. |
| `graficar_diversidad` | `graficar_diversidad(historial; titulo="Apellidos vivos por generacion")` | Nº de apellidos vivos por generación. |
| `graficar_dist_n_vivos` | `graficar_dist_n_vivos(rep; titulo=...)` | Barras: fracción de corridas que terminaron con k apellidos vivos. |
| `graficar_prob_ganar` | `graficar_prob_ganar(rep; titulo=...)` | Barras: probabilidad de que cada apellido sea el único superviviente. |

Función privada: `_mapa_indices(compatibles)` (mapea ids de apellido a índices
de vértice para Graphs.jl).

## Dependencias

- Externas: `Plots`, `Graphs`, `GraphRecipes`, `StatsBase`.
- Del proyecto: `src/core/experimentMetrics.jl` (usa `conteos_por_gen` y
  `apellidos_vivos_por_gen`), que debe estar cargado antes.

## Uso

```julia
include(joinpath("src", "coupleModel.jl"))
include(joinpath("src", "core", "graphTopologies.jl"))
include(joinpath("src", "core", "experimentMetrics.jl"))
include(joinpath("src", "core", "visualization.jl"))

compat = topo_cadena(4)
res = simular_parejas(poblacion_desde_parejas(Dict(a => 3 for a in 1:4)), compat)
dibujar_grafo(compat)
graficar_diversidad(res.historial)
```
