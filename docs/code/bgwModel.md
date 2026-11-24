# `src/bgwModel.jl` — Modelo 1: Bienaymé-Galton-Watson

## Qué hace

Modela un apellido como un linaje bajo un proceso de ramificación de
Bienaymé-Galton-Watson (BGW). Cada portador tiene un número aleatorio de hijos
(que heredan el apellido) según una distribución de descendencia `pk`. Calcula
la probabilidad de extinción `q` del linaje y simula su evolución generación a
generación.

## Cómo lo hace

Sea `Z_n` el número de portadores en la generación `n`. Entonces
`Z_{n+1} = X_1 + ... + X_{Z_n}`, con cada `X_i` i.i.d. de distribución `pk`. Si
`Z_n = 0`, extinción (estado absorbente).

La probabilidad de extinción `q` es el menor punto fijo de la función generadora
de probabilidad `g(t) = sum(pk[k] t^(k-1))`, es decir la menor raíz de
`g(t) = t`. Se obtiene buscando la raíz de `h(t) = g(t) - t` en `[0, 1-margen]`
con los métodos de `src/core/numericalMethods.jl` (evitando la raíz trivial
`t = 1`). El resultado se valida empíricamente por Monte Carlo.

Regímenes según la media `m = E[X]`: subcrítico (`m < 1`) y crítico (`m = 1`)
dan extinción segura (`q = 1`); supercrítico (`m > 1`) permite supervivencia
(`q < 1`).

## Qué expone

| Símbolo | Firma | Propósito |
|:--|:--|:--|
| `ModeloPBGW` | `struct` con campo `pk::Vector{Float64}` | Distribución de descendencia. |
| `ModeloPBGW` | `ModeloPBGW(pk::Vector{Float64})` | Construye desde un vector `pk` (valida no negatividad y suma 1). |
| `ModeloPBGW` | `ModeloPBGW(tasa_natalidad::Float64; max_hijos=15)` | Construye desde una Poisson de media `tasa_natalidad`. |
| `esperanza` | `esperanza(modelo)` | Media de descendencia `m = E[X]`. |
| `es_extincion_segura` | `es_extincion_segura(modelo)` | `true` si `m <= 1`. |
| `fgp` | `fgp(modelo, t::Float64)` | Evalúa la función generadora `g(t)`. |
| `probabilidad_extincion` | `probabilidad_extincion(modelo; buscador=:biseccion)` | Calcula `q` (buscador `:biseccion` o `:secante`). |
| `simular_generaciones` | `simular_generaciones(modelo, generaciones; poblacion_inicial=1, tope_poblacion=1_000_000, rng=...)` | Trayectoria estocástica de un linaje (`Vector{Int}` de longitud `generaciones+1`). |
| `frecuencia_extincion` | `frecuencia_extincion(modelo, generaciones; n_sim=10_000, poblacion_inicial=1, tope_poblacion=1_000_000, rng=...)` | Estimación Monte Carlo de `q`. |

`tope_poblacion` acota los linajes supercríticos que crecen exponencialmente: al
superarlo se consideran "establecidos" (ya no se extinguirán) y la población se
fija al tope, para que la simulación siga siendo viable.

## Dependencias

- Externas: `Distributions`, `Random`.
- Del proyecto: `src/core/numericalMethods.jl` (usa `biseccion` y `secante`),
  que debe estar cargado antes.

## Uso

```julia
include(joinpath("src", "core", "numericalMethods.jl"))
include(joinpath("src", "bgwModel.jl"))

m = ModeloPBGW(2.0)                    # Poisson(λ=2)  → supercrítico
q = probabilidad_extincion(m)          # probabilidad de extinción teórica
traj = simular_generaciones(m, 50)     # una trayectoria de 50 generaciones
f = frecuencia_extincion(m, 50)        # q estimada por Monte Carlo
```
