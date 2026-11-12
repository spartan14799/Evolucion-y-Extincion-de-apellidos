# Modelo estacionario (población de tamaño constante)

Este documento describe el **modelo estacionario** de evolución y extinción de
apellidos implementado en [`src/stacional-model.jl`](../src/stacional-model.jl)
y ejecutado desde [`scripts/run_stacional_model.jl`](../scripts/run_stacional_model.jl).

A diferencia del proceso de Bienaymé–Galton–Watson (ver
[`docs/justificacion-teorica.md`](justificacion-teorica.md)), donde la población
puede crecer o decrecer, aquí el tamaño de la población **se conserva en cada
generación**.

## 1. Idea del modelo

- La población tiene un tamaño fijo `N` que **no cambia** entre generaciones.
- En cada generación los individuos se **mezclan al azar** y se agrupan en
  parejas (bloques) de tamaño `k` (por defecto `k = 2`).
- De cada grupo **sobrevive un único apellido** (el del primer miembro tras la
  mezcla), que se propaga a todos los miembros del grupo en la siguiente
  generación.
- Como consecuencia, la **diversidad de apellidos nunca aumenta**: solo puede
  mantenerse o disminuir. El proceso avanza siempre hacia la **fijación** de un
  único apellido, que es un **estado absorbente**.

Este comportamiento es análogo a la **deriva genética**: por puro azar del
muestreo, una población finita termina perdiendo variantes hasta quedarse con
una sola.

> Nota: para que el tamaño se conserve, `N` debe ser **divisible por `k`**. Con
> `k = 2` esto significa que la población debe ser **par**. El código valida
> esto y lanza un error claro si no se cumple.

## 2. Mapeo al código

| Concepto del modelo                 | Elemento en el código                          |
| ----------------------------------- | ---------------------------------------------- |
| Estado de la población              | `Vector{Int}` de ids de apellido               |
| Parámetros de la dinámica           | `struct ParametrosModelo`                      |
| Una transición de generación        | `siguiente_generacion(poblacion, rng; parametros)` |
| Simulación hasta fijación/límite    | `simular_fijacion(estado_inicial; ...)`        |
| Construir población inicial         | `configurar_poblacion(...)` (en el script)     |

### `ParametrosModelo`

Agrupa los parámetros configurables de la dinámica:

| Campo              | Tipo  | Por defecto | Significado                                             |
| ------------------ | ----- | ----------- | ------------------------------------------------------- |
| `tamano_pareja`    | `Int` | `2`         | Individuos por grupo que eligen un único apellido.      |
| `max_generaciones` | `Int` | `10_000`    | Tope de seguridad de generaciones antes de rendirse.    |

```julia
# Parámetros por defecto
p = ParametrosModelo()

# Parejas de 2, máximo 5000 generaciones
p = ParametrosModelo(tamano_pareja = 2, max_generaciones = 5_000)
```

### `siguiente_generacion(poblacion, rng; parametros)`

Aplica **una** transición:

1. Mezcla la población al azar con `shuffle`.
2. La recorre en bloques de `tamano_pareja`.
3. Copia el apellido del primer miembro de cada bloque a todo el bloque.

Devuelve un nuevo vector del **mismo tamaño**. Lanza error si el tamaño no es
divisible por `tamano_pareja`.

### `simular_fijacion(estado_inicial; parametros, max_generaciones, rng)`

Itera `siguiente_generacion` hasta que quede **un único apellido** o se alcance
`max_generaciones`. Devuelve una *named tuple*:

| Campo          | Significado                                                       |
| -------------- | ---------------------------------------------------------------- |
| `ganador`      | Id del apellido fijado, o `nothing` si no hubo fijación.         |
| `generaciones` | Número de transiciones simuladas.                                |
| `historial`    | Vector con el estado completo de cada generación (incluye `0`).  |
| `fijacion`     | `true` si se alcanzó un único apellido.                          |

Invariantes garantizadas (verificadas en los tests):

- `length(historial) == generaciones + 1`.
- Todos los estados del `historial` tienen el mismo tamaño que `estado_inicial`.
- Si `fijacion == true`, todo el estado final es igual a `ganador`.

## 3. Cómo cambiar los parámetros

Hay dos formas de configurar una ejecución.

### a) Editando el script

Abre [`scripts/run_stacional_model.jl`](../scripts/run_stacional_model.jl) y
modifica el bloque `PARÁMETROS DE LA SIMULACIÓN` dentro de `main()`:

```julia
POBLACION_INICIAL   = 10      # Número de individuos (debe ser divisible por TAMANO_PAREJA)
APELLIDOS_INICIALES = 4       # Apellidos distintos iniciales (<= población)
TAMANO_PAREJA       = 2       # Individuos por grupo
MAX_GENERACIONES    = 10_000  # Tope de seguridad
SEMILLA             = 42      # Entero fijo (reproducible) o `nothing` (azar)
```

### b) Con variables de entorno (sin tocar el código)

Las variables de entorno **tienen prioridad** sobre los valores del código:

```bash
POBLACION_INICIAL=20 APELLIDOS_INICIALES=6 SEMILLA=7 \
    julia --project=. scripts/run_stacional_model.jl
```

Variables reconocidas:

| Variable              | Por defecto | Significado                                          |
| --------------------- | ----------- | ---------------------------------------------------- |
| `POBLACION_INICIAL`   | `10`        | Número de individuos (divisible por `TAMANO_PAREJA`).|
| `APELLIDOS_INICIALES` | `4`         | Apellidos distintos iniciales (`<=` población).      |
| `TAMANO_PAREJA`       | `2`         | Individuos por grupo que eligen un apellido.         |
| `MAX_GENERACIONES`    | `10000`     | Tope de generaciones de seguridad.                   |
| `SEMILLA`             | `42`        | Entero para reproducibilidad, o `nothing` para azar. |

### c) Desde tu propio código Julia

```julia
include("src/stacional-model.jl")
using Random

parametros = ParametrosModelo(tamano_pareja = 2, max_generaciones = 2_000)
estado = [1, 1, 2, 3, 4, 4, 5, 6]          # 8 individuos, 6 apellidos distintos
res = simular_fijacion(estado; parametros = parametros, rng = MersenneTwister(42))

println(res.fijacion)       # true si se fijó un apellido
println(res.ganador)        # id del apellido ganador
println(res.generaciones)   # cuántas generaciones tardó
```

## 4. Ejecutar la simulación

Desde la raíz del proyecto:

```bash
julia --project=. scripts/run_stacional_model.jl
```

La salida muestra el estado y las frecuencias de apellidos en cada generación, y
termina indicando si se alcanzó la fijación (y con qué apellido) o si se agotó el
límite de generaciones.

## 5. Ejecutar los tests

```bash
julia --project=. test/test_stacional.jl
```

Los tests (en [`test/test_stacional.jl`](../test/test_stacional.jl)) cubren:

1. **`ParametrosModelo`**: valores por defecto y constructor por palabras clave.
2. **Conservación del tamaño**: `siguiente_generacion` nunca cambia `N` y la
   diversidad no aumenta.
3. **Tamaño de pareja configurable**: funciona con `k = 4` y rechaza tamaños no
   divisibles o parejas inválidas.
4. **Reproducibilidad**: misma semilla ⇒ mismo resultado.
5. **Estado absorbente y límites**: población ya fijada, estructura del
   historial y comportamiento con `max_generaciones = 0`.
6. **Convergencia a fijación**: para poblaciones pequeñas, con un tope amplio la
   fijación se alcanza siempre (probado con 20 semillas distintas).

## 6. Notas y limitaciones

- **La población debe ser divisible por `tamano_pareja`.** El modelo original
  no protegía este caso: con una población impar y parejas de 2, el tamaño de la
  población *cambiaba* entre generaciones, rompiendo la hipótesis de
  estacionariedad. Ahora `siguiente_generacion` lanza un error explícito.
- La elección del "ganador" es el **primer** elemento de cada bloque tras la
  mezcla aleatoria; dado que la mezcla es uniforme, cada apellido del bloque
  tiene la misma probabilidad de ganar.
- El modelo no incorpora mutación ni aparición de nuevos apellidos: la
  diversidad solo puede disminuir.
