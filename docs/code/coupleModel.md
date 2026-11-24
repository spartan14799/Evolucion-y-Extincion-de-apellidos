# `src/coupleModel.jl` — Modelo 2: modelo por parejas

## Qué hace

Modela la evolución de apellidos sobre un grafo donde cada nodo es un apellido
y las aristas (`compatibles`) indican qué apellidos pueden emparejarse. La
unidad reproductiva es la **pareja**: cada pareja aporta exactamente 2 hijos, de
modo que el tamaño total `N = 2P` se conserva. Simula la evolución hasta
fijación, equilibrio permanente o un tope de generaciones. Es un motor
independiente (no usa otros archivos del proyecto).

## Cómo lo hace

Una generación transforma `2P` padres en `2P` hijos:

1. **Emparejamiento sin reemplazo.** Se barajan las personas y, en orden
   aleatorio, cada persona suelta busca pareja entre las personas todavía
   sueltas cuyo apellido sea compatible. La pareja concreta se sortea según el
   modo de probabilidad, con pesos calculados una vez al inicio de la generación.
2. **Herencia.** Por cada pareja se elige al azar (50/50) a uno de los dos
   miembros; ambos hijos toman ese apellido.
3. Cada pareja produce 2 hijos, así `N = 2P` se conserva.

Si al final quedan personas sueltas sin ningún compatible disponible, se
emparejan de forma **forzada** (incompatibles); cada evento es una *mutación
definida* y se registra.

La diversidad sólo baja cuando se forma una pareja **mixta** (dos apellidos
distintos compatibles). Si entre los apellidos vivos no existe ninguna arista
`i -> j` con `i != j`, ya no hay parejas mixtas posibles: el estado es un
equilibrio permanente (criterio exacto `hay_parejas_mixtas_posibles`). La
simulación para al fijar un apellido, al llegar a ese equilibrio, o al tope.

Modos de probabilidad (`modo_prob`):
- `:uniforme_aristas` (default): cada persona elige destino uniforme entre sus
  aristas hacia apellidos vivos.
- `:ponderado`: probabilidad proporcional al conteo del apellido destino.

## Qué expone

| Símbolo | Firma | Propósito |
|:--|:--|:--|
| `ParametrosParejas` | `struct` (kwdef) `max_generaciones::Int=1000`, `modo_prob::Symbol=:uniforme_aristas` | Parámetros de la simulación. |
| `ResultadoParejas` | `struct` con `estado_final`, `generaciones`, `historial`, `forzadas_por_gen`, `absorbido`, `apellidos_final`, `transitorio` | Resultado de una corrida. |
| `contar` | `contar(poblacion::Vector{Int}) -> Dict{Int,Int}` | Conteo de personas por apellido. |
| `pesos_emparejamiento` | `pesos_emparejamiento(apellido_i, compatibles, conteos, vivos, modo) -> (destinos, pesos)` | Destinos posibles y sus pesos para una persona. |
| `siguiente_generacion_parejas` | `siguiente_generacion_parejas(poblacion, compatibles, rng=...; modo_prob=:uniforme_aristas) -> (nueva_poblacion, parejas_forzadas)` | Una transición generacional. |
| `hay_parejas_mixtas_posibles` | `hay_parejas_mixtas_posibles(apellidos_vivos, compatibles) -> Bool` | Criterio exacto de progreso de la diversidad. |
| `simular_parejas` | `simular_parejas(estado_inicial, compatibles; parametros=ParametrosParejas(), rng=...) -> ResultadoParejas` | Corrida completa hasta absorción o tope. |

Funciones privadas (prefijo `_`): `_elegir_companero!`, `_tomar_cualquier_suelto!`.

Sobre `ResultadoParejas`: `absorbido=true` indica fijación de un apellido o
equilibrio permanente; `transitorio=true` indica que se cortó por
`max_generaciones` y la diversidad aún podría bajar (excluyente con `absorbido`).

## Dependencias

- Externas: `Random`, `StatsBase`.
- Del proyecto: ninguna (motor independiente). Las topologías de `compatibles`
  suelen venir de `src/core/graphTopologies.jl`.

## Uso

```julia
include(joinpath("src", "coupleModel.jl"))

estado = [1, 1, 2, 2]                      # 2 apellidos, 2 personas cada uno
compatibles = Dict(1 => [1, 2], 2 => [1, 2])
res = simular_parejas(estado, compatibles; parametros = ParametrosParejas(max_generaciones = 500))
res.absorbido, res.apellidos_final        # ¿absorbió?, apellidos vivos al final
```
