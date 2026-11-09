# Justificación teórica: evolución y extinción de apellidos

Este documento describe el modelo matemático que sustenta la simulación, cómo
se traduce en código y por qué los métodos numéricos del proyecto son la pieza
central del cálculo.

---

## El problema

La transmisión de un apellido sigue, de forma simplificada, la misma lógica que
la transmisión de un linaje por línea paterna: cada portador del apellido tiene
un número aleatorio de descendientes que heredan ese apellido. Generación tras
generación, algunos linajes crecen y otros desaparecen.

La pregunta central es:

> **¿Cuál es la probabilidad de que un apellido termine extinguiéndose?**

Para responderla se usa un modelo clásico de probabilidad: el **proceso de
Bienaymé–Galton–Watson (BGW)**, también conocido como proceso de ramificación.

---

## El proceso de Galton-Watson

### Definición

Consideramos una población de portadores de un apellido que evoluciona en
generaciones discretas `n = 0, 1, 2, ...`. Sea `Z_n` el número de portadores en
la generación `n`. El proceso se define así:

- Se parte de `Z_0` portadores (habitualmente `Z_0 = 1`).
- Cada portador `i` de la generación `n` produce un número aleatorio de hijos
  `X_i`, independiente e idénticamente distribuido según una ley común.
- La siguiente generación es la suma de todos los hijos:

$$Z_{n+1} = \sum_{i=1}^{Z_n} X_i$$

- Si `Z_n = 0`, el apellido se ha **extinguido** y permanece en 0 para siempre
  (estado absorbente).

### La distribución de descendencia

La ley común de `X_i` se describe mediante la secuencia de probabilidades:

$$p_k = P(X = k), \quad k = 0, 1, 2, \dots$$

donde `p_k` es la probabilidad de que un portador tenga exactamente `k` hijos.
Se cumple que todas son no negativas y suman 1.

En el código, esto es el vector `pk` de la estructura `ModeloPBGW`
(`src/markovModels.jl`). Por convención de indexado en Julia:

- `pk[1]` = P(0 hijos)
- `pk[2]` = P(1 hijo)
- `pk[k]` = P(k−1 hijos)

El modelo admite dos formas de construcción:

1. **Directa**: se pasa el vector `pk` explícito.
2. **Por tasa de natalidad**: se pasa `λ` y se usa una distribución de **Poisson**
   truncada y normalizada, `X ~ Poisson(λ)`. Poisson es una elección natural para
   modelar el número de descendientes.

---

## La esperanza y los tres regímenes

El parámetro más importante del proceso es la **media de descendencia**:

$$m = E[X] = \sum_{k=0}^{\infty} k \, p_k$$

En el código es la función `esperanza(modelo)`. Este número determina por
completo el destino cualitativo del linaje:

| Régimen | Condición | Comportamiento | Prob. de extinción `q` |
|---|---|---|---|
| **Subcrítico** | `m < 1` | Se extingue casi seguro | `q = 1` |
| **Crítico** | `m = 1` | Se extingue casi seguro* | `q = 1` |
| **Supercrítico** | `m > 1` | Puede sobrevivir | `q < 1` |

\* En el caso crítico la extinción también es segura (salvo el caso degenerado
`p_1 = 1`), aunque el tiempo medio hasta la extinción es infinito.

La función `es_extincion_segura(modelo)` devuelve `true` cuando `m ≤ 1`.

Interpretación intuitiva: si en promedio cada persona tiene menos de un hijo que
herede el apellido, el linaje se contrae y acaba desapareciendo. Solo cuando el
promedio supera 1 el apellido tiene una probabilidad positiva de perdurar
indefinidamente.

---

## Cálculo de la probabilidad de extinción

### 4.1 La función generadora de probabilidad

La herramienta clave es la **función generadora de probabilidad** (FGP) de la
distribución de descendencia:

$$g(t) = E[t^X] = \sum_{k=0}^{\infty} p_k \, t^k, \quad t \in [0, 1]$$

En el código es `fgp(modelo, t)`.

Propiedades que usaremos:

- `g(0) = p_0` (la probabilidad de no tener hijos).
- `g(1) = 1` (las probabilidades suman 1).
- `g` es continua, creciente y convexa en `[0, 1]`.
- `g'(1) = m` (la derivada en 1 es la media de descendencia).

### La ecuación de punto fijo

El resultado fundamental de la teoría de procesos de ramificación establece que
la **probabilidad de extinción `q` es el menor punto fijo de `g` en `[0, 1]`**,
es decir, la menor solución de:

$$g(q) = q$$

- `t = 1` siempre es solución (porque `g(1) = 1`); es la **raíz trivial**.
- Si `m > 1`, existe además una única solución `q` en el intervalo abierto
  `(0, 1)`, y **esa** es la probabilidad de extinción.
- Si `m ≤ 1`, la única solución en `[0, 1]` es `q = 1`.

### Conexión con los métodos numéricos

Encontrar `q` equivale a **resolver una ecuación no lineal**, que es exactamente
el propósito de `src/numericalMethods.jl`. Reescribimos el punto fijo como una
búsqueda de raíz definiendo:

$$h(t) = g(t) - t$$

Entonces `q` es la raíz de `h(t) = 0` en `(0, 1)`. La función
`probabilidad_extincion(modelo)` procede así:

1. Si `m ≤ 1`, devuelve `q = 1` directamente (no hace falta resolver nada).
2. Si `m > 1`, construye `h(t) = fgp(modelo, t) - t` y busca su raíz en
   `[0, 1 − ε]` usando los métodos numéricos del proyecto:
   - **Bisección** (`biseccion`) por defecto: robusta, confinada al intervalo.
   - **Secante** (`secante`) como opción alternativa.

> **Nota sobre la secante.** La secante es un método *abierto*: no queda
> confinado al intervalo de búsqueda. Como `t = 1` es siempre raíz de `h`,
> arrancar cerca de 1 puede hacer que la secante converja a esa raíz trivial en
> lugar de a `q`. Por eso `probabilidad_extincion` primero acota `q` con unos
> pasos de bisección y luego lanza la secante en un intervalo estrecho alrededor
> de la aproximación, donde el único punto fijo es `q`.

### Ejemplo con solución analítica

Para validar, consideremos `p_0 = 0.4`, `p_1 = 0`, `p_2 = 0.6` (es decir,
`pk = [0.4, 0.0, 0.6]`). Aquí `m = 2 · 0.6 = 1.2 > 1`, luego es supercrítico.

La ecuación de punto fijo es:

$$g(t) = 0.4 + 0.6\,t^2 = t \;\Longrightarrow\; 0.6\,t^2 - t + 0.4 = 0$$

Resolviendo la cuadrática:

$$t = \frac{1 \pm \sqrt{1 - 0.96}}{1.2} = \frac{1 \pm 0.2}{1.2}$$

Las dos raíces son `t = 1` (trivial) y `t = 0.8/1.2 = 2/3 ≈ 0.6667`. La menor es
`q = 2/3`, que es la probabilidad de extinción. El test
`test/test_markov.jl` comprueba que tanto bisección como secante recuperan este
valor.

---

## Simulación estocástica (Monte Carlo)

Además del cálculo teórico, el proyecto incluye una **simulación directa** del
proceso, que sirve como validación independiente.

### Una trayectoria

`simular_generaciones(modelo, generaciones)` simula un linaje concreto:

- Parte de `poblacion_inicial` portadores.
- En cada generación, muestrea el número de hijos de cada portador a partir de
  `pk` (usando una distribución `Categorical`) y suma los resultados.
- Devuelve el vector de población por generación. Si llega a 0, el linaje muere.

**Tope de población.** En modelos supercríticos, los linajes que *no* se
extinguen crecen exponencialmente (a ritmo `m` por generación), lo que haría la
simulación inviable (millones de individuos tras unas pocas decenas de
generaciones). Para evitarlo se usa el parámetro `tope_poblacion`: cuando un
linaje supera ese umbral se considera **establecido** (ya no se extinguirá) y la
población se fija al tope. Esto no altera el resultado de extinción y mantiene el
coste acotado.

### Estimación de la probabilidad de extinción

`frecuencia_extincion(modelo, generaciones; n_sim)` ejecuta `n_sim` simulaciones
independientes y devuelve la fracción de linajes que terminaron extintos. Es una
estimación empírica de `q` por el método de Monte Carlo.

### Validación cruzada

El diseño permite obtener `q` por dos caminos distintos:

1. **Teórico / numérico**: resolver `g(t) = t` con bisección → valor exacto.
2. **Empírico / simulado**: correr el proceso muchas veces y contar extinciones.

Si ambos coinciden dentro del margen estadístico, es una fuerte evidencia de que
tanto el modelo como los métodos numéricos están correctamente implementados.
Resultados observados en el proyecto:

| Modelo | `q` teórica (bisección) | `q` simulada (Monte Carlo) |
|---|---|---|
| `pk = [0.4, 0, 0.6]` (m = 1.2) | 0.6666 | ≈ 0.6647 |
| Poisson(λ = 1.5) (m = 1.5) | 0.4172 | ≈ 0.4128 |
| Poisson(λ = 0.8) (m = 0.8) | 1.0000 | 1.0000 |

---

## Mapa del código

| Concepto matemático | Función / tipo | Archivo |
|---|---|---|
| Distribución de descendencia `p_k` | `ModeloPBGW` | `src/markovModels.jl` |
| Media `m = E[X]` | `esperanza` | `src/markovModels.jl` |
| Criterio `m ≤ 1` | `es_extincion_segura` | `src/markovModels.jl` |
| FGP `g(t)` | `fgp` | `src/markovModels.jl` |
| Probabilidad de extinción `q` | `probabilidad_extincion` | `src/markovModels.jl` |
| Simulación de una trayectoria | `simular_generaciones` | `src/markovModels.jl` |
| Estimación Monte Carlo de `q` | `frecuencia_extincion` | `src/markovModels.jl` |
| Búsqueda de raíces | `biseccion`, `newton`, `secante` | `src/numericalMethods.jl` |

El flujo esencial: `probabilidad_extincion` envuelve `fgp` en `h(t) = g(t) − t`
y se la entrega a `biseccion` / `secante` de `numericalMethods.jl` para hallar
la raíz, que es la probabilidad de extinción buscada.

- La distribución de Poisson como modelo del número de descendientes.
