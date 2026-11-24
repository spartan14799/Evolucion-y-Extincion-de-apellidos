# `test/` — Tests y benchmark

Tres archivos ejecutables desde la raíz con `julia --project=. test/...` (o vía
`make`). Cubren los métodos numéricos y el modelo BGW.

> Nota: el modelo por parejas (`coupleModel`, `experimentMetrics`,
> `graphTopologies`) no tiene tests automatizados; se valida ejecutando
> `scripts/experimentos_parejas.jl` y los notebooks.

---

## `test/runtests.jl` — Tests de los métodos numéricos

**Qué prueba.** La correctitud de `biseccion`, `newton` y `secante` sobre
funciones con raíz conocida (p. ej. `x²-2`, `cos(x)-x`, `x³-x-2`): verifica que
la raíz encontrada esté dentro de tolerancia, que `convergio` sea `true`, la
consistencia entre los tres métodos, y los casos de error (intervalo sin cambio
de signo en bisección, derivada nula en Newton).

**Dependencias.** `Test` (stdlib). Incluye `src/core/numericalMethods.jl`.

**Uso.** `julia --project=. test/runtests.jl` (o `make test`). Esperado: 25/25.

---

## `test/test_markov.jl` — Tests del modelo BGW

**Qué prueba.** El modelo de Bienaymé-Galton-Watson en 7 bloques con progreso en
vivo: construcción de `ModeloPBGW` (y validaciones del constructor), `esperanza`,
`es_extincion_segura`, `fgp`, `probabilidad_extincion` (bisección y secante, con
un caso analítico `q = 2/3`), `simular_generaciones`, y una validación Monte
Carlo de `frecuencia_extincion` frente a la `q` teórica.

**Dependencias.** `Test`, `Random`, `Printf`. Incluye
`src/core/numericalMethods.jl` y `src/bgwModel.jl`.

**Uso.** `julia --project=. test/test_markov.jl` (o `make test-markov`).
Esperado: 28/28. El bloque de Monte Carlo tarda unos segundos.

---

## `test/benchmarks.jl` — Benchmark de los métodos numéricos

**Qué hace.** Mide el rendimiento de `biseccion`, `newton` y `secante`. Usa
`BenchmarkTools` si está instalado (estadísticas robustas); si no, recurre a
`@elapsed` de la librería estándar, para funcionar sin dependencias.

**Dependencias.** `Printf`; `BenchmarkTools` opcional. Incluye
`src/core/numericalMethods.jl`.

**Uso.** `julia --project=. test/benchmarks.jl` (o `make bench`).
