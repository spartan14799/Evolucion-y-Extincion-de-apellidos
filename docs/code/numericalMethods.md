# `src/core/numericalMethods.jl` — Métodos numéricos de búsqueda de raíces

## Qué hace

Implementa tres métodos clásicos de búsqueda de raíces de una función real:
bisección, Newton y secante. Es soporte genérico; el modelo BGW los usa para
resolver `g(t) = t` y hallar la probabilidad de extinción.

## Cómo lo hace

- **Bisección:** método cerrado. Requiere que `f` tenga signos opuestos en los
  extremos `[a, b]` y reduce el intervalo a la mitad en cada paso. Garantiza
  convergencia si hay cambio de signo.
- **Newton:** método abierto. Usa `f` y su derivada `df` desde un punto inicial
  `x0`: `x_{k+1} = x_k - f(x_k)/f'(x_k)`. Falla si la derivada es cero.
- **Secante:** método abierto. Aproxima la derivada con dos puntos; reordena
  para que `a` sea la mejor aproximación en cada iteración.

Las tres devuelven un `NamedTuple` `(raiz, iteraciones, convergio)` y emiten un
`@warn` si no convergen dentro de `M` iteraciones.

## Qué expone

| Símbolo | Firma | Propósito |
|:--|:--|:--|
| `biseccion` | `biseccion(f::Function, a::Float64, b::Float64; M=100, delta=1e-5, epsilon=1e-5)` | Raíz por bisección en `[a, b]`. |
| `newton` | `newton(f::Function, df::Function, x0::Float64; M=50, delta=1e-5, epsilon=1e-5)` | Raíz por Newton desde `x0` usando la derivada `df`. |
| `secante` | `secante(f::Function, a::Float64, b::Float64; M=50, delta=1e-5, epsilon=1e-5)` | Raíz por secante desde `a`, `b`. |

Criterios de parada: `delta` (cambio en la variable) y `epsilon` (valor de `f`).

## Dependencias

Ninguna (librería estándar de Julia).

## Uso

```julia
include(joinpath("src", "core", "numericalMethods.jl"))

f(x) = x^2 - 2
r = biseccion(f, 0.0, 2.0)        # r.raiz ≈ 1.41421
rn = newton(f, x -> 2x, 1.0)      # con derivada explícita
rs = secante(f, 0.0, 2.0)
```
