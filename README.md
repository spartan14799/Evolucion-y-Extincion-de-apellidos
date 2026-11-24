# Evolución y extinción de apellidos

Simulación y análisis de la **probabilidad de extinción de un apellido** usando
un proceso de Bienaymé–Galton–Watson, resuelto con métodos numéricos clásicos
(bisección, Newton, secante) y validado mediante simulación Monte Carlo.

El proyecto está escrito en [Julia](https://julialang.org/).

## ¿Qué hace?

- Modela la descendencia de un apellido como un proceso de ramificación.
- Calcula la probabilidad de extinción `q` resolviendo la ecuación de punto fijo
  `g(t) = t` con métodos numéricos de búsqueda de raíces.
- Simula la evolución de linajes generación a generación.
- Contrasta el resultado teórico contra una estimación Monte Carlo.

La teoría del proceso de ramificación está en
[`docs/justificacion-teorica.md`](docs/justificacion-teorica.md). El modelo por
parejas sobre el grafo de apellidos (reglas, modos de probabilidad, criterio de
parada y catálogo de experimentos) está en
[`docs/modelo-parejas.md`](docs/modelo-parejas.md).

## Estructura

```
.
├── Project.toml              # Dependencias del proyecto (entorno aislado)
├── Makefile                  # Atajos de ejecución (make help)
├── README.md
├── docs/
│   ├── justificacion-teorica.md   # Modelo matemático y su mapeo al código
│   ├── modelo-parejas.md          # Modelo por parejas: reglas, parada y experimentos
│   └── code/                      # Una ficha por archivo de código (ver docs/code/README.md)
├── src/
│   ├── bgwModel.jl           # MODELO 1: Bienaymé–Galton–Watson (extinción de linajes)
│   ├── coupleModel.jl        # MODELO 2: modelo por parejas sobre el grafo de apellidos
│   ├── core/                 # Requisitos genéricos de los modelos
│   │   ├── numericalMethods.jl   # biseccion, newton, secante
│   │   ├── graphTopologies.jl    # Generadores de topologías del grafo de apellidos
│   │   ├── experimentMetrics.jl  # Métricas y runner de experimentos del modelo por parejas
│   │   └── visualization.jl      # Utilidades de graficación (grafo, conteos, diversidad)
│   └── experiments/          # Datos reales y orquestación para los notebooks de results/
│       ├── daneData.jl           # Datos del Censo DANE 2018 (solo datos)
│       ├── daneBgwRunners.jl     # Monte Carlo BGW sobre los datos DANE
│       └── daneCoupleRunners.jl  # Árbol genealógico por parejas sobre los datos DANE
├── scripts/
│   ├── pluto.jl                  # Lanzador multiplataforma de los notebooks Pluto
│   ├── run_simulation.jl         # Simulación inicial de demostración (BGW)
│   └── experimentos_parejas.jl   # Catálogo reproducible: escribe PNG + CSV en output/
├── notebooks/
│   ├── basic/                # Ilustración interactiva de cada modelo
│   │   ├── bgw_model.jl          # Modelo BGW: regímenes, FGP, Monte Carlo, cadena de Markov
│   │   └── couple_model.jl       # Modelo por parejas: topologías y dinámica del grafo
│   └── results/              # Aplicación a datos reales (Censo DANE 2018)
│       ├── bgw_colombia_dane.jl       # ¿Cuándo se extinguen los apellidos de Colombia?
│       └── couple_genealogy_dane.jl   # Árbol genealógico por parejas con datos DANE
└── test/
    ├── runtests.jl           # Tests de los métodos numéricos
    ├── test_markov.jl        # Tests del modelo BGW / simulación
    └── benchmarks.jl         # Benchmark de los métodos numéricos
```

La documentación detallada de cada archivo de código (qué hace, cómo, qué
funciones expone y su uso) está en [`docs/code/`](docs/code/README.md). El
encabezado de cada archivo `.jl` lleva sólo una descripción breve y apunta a su
ficha correspondiente.

## Requisitos

- Julia **1.12**.

Comprueba tu versión con:

```bash
julia --version
```

### Versión de Julia fija en todos los dispositivos

El proyecto se desarrolla y prueba con **Julia 1.12**, y el `Manifest.toml`
versionado en git está resuelto con esa versión. Usar una versión distinta
(p. ej. 1.13) provoca errores de precompilación porque el `Manifest.toml`
deja de ser compatible.

Para garantizar la misma versión en cualquier máquina usamos
[`juliaup`](https://github.com/JuliaLang/juliaup) (el gestor oficial de
versiones de Julia) con un *override* de carpeta: así, dentro de este
directorio, `julia` siempre apunta a 1.12 aunque tu versión por defecto sea
otra.

En **cada dispositivo nuevo**, tras clonar el repo, ejecuta una vez:

```bash
juliaup add 1.12            # instala Julia 1.12 si no la tienes
juliaup override set 1.12   # fija 1.12 solo para esta carpeta
```

El *override* se guarda en la configuración de juliaup de cada máquina (no se
versiona en git), por eso hay que repetir estos dos comandos una vez por
dispositivo. Verifica el estado con:

```bash
juliaup override status     # lista los overrides por carpeta
julia --version             # debe decir 1.12.x dentro de esta carpeta
```

> Si no usas juliaup, basta con tener instalado Julia 1.12 y asegurarte de que
> es el `julia` de tu `PATH` al trabajar en el proyecto.

### Pipeline de puesta a punto (resumen)

```bash
git clone https://github.com/spartan14799/Evolucion-y-Extincion-de-apellidos.git && cd Evolucion-y-Extincion-de-apellidos
juliaup add 1.12 && juliaup override set 1.12   # misma versión de Julia
make instantiate                                 # instala las deps EXACTAS del Manifest.toml
make pluto                                       # lanza los notebooks
```

`make instantiate` lee el `Manifest.toml` versionado e instala exactamente esas
versiones de paquetes, por lo que el entorno queda idéntico en todos los
dispositivos siempre que la versión de Julia coincida (de ahí el *override*).

## Atajos con `make` para la ejecución

El proyecto incluye un `Makefile` que envuelve todos los comandos habituales.
Ejecuta `make` (o `make help`) para ver la lista:

```bash
make help
```

Objetivos principales:

| Comando            | Qué hace                                                  |
| ------------------ | --------------------------------------------------------- |
| `make instantiate` | Instala y precompila las dependencias (solo una vez)      |
| `make pluto`       | Lanza el servidor Pluto (eliges el notebook en el navegador) |
| `make bgw`         | Abre el notebook básico del modelo Bienaymé–Galton–Watson |
| `make parejas`     | Abre el notebook básico del modelo por parejas            |
| `make bgw-dane`    | Abre el notebook de resultados BGW con datos del Censo DANE 2018 |
| `make parejas-dane`| Abre el notebook de resultados del árbol genealógico (DANE 2018) |
| `make experimentos`| Corre el catálogo de experimentos por parejas (PNG + CSV en `output/`) |
| `make simulation`  | Ejecuta la simulación de demostración                     |
| `make test`        | Tests de los métodos numéricos                            |
| `make test-markov` | Tests del modelo de Markov / simulación                   |
| `make bench`       | Benchmark de los métodos numéricos                        |

El `Makefile` usa el `julia` de tu `PATH`. Si lo tienes instalado en otra ruta
o con otro nombre, pásalo como variable: `make pluto JULIA=/ruta/a/julia`.

> **Windows:** `make` no viene preinstalado. Si no lo tienes, usa directamente
> los comandos `julia` equivalentes (cada fila de arriba corresponde al comando
> `julia --project=. ...` documentado en las secciones siguientes). Los antiguos
> archivos `.bat` se eliminaron porque ataban la ruta de Julia a una versión
> concreta; el lanzador `scripts/pluto.jl` los reemplaza de forma portable.

## Notebooks de Pluto

Los notebooks interactivos viven en `notebooks/` y se abren con Pluto. El
lanzador `scripts/pluto.jl` es multiplataforma: usa el `julia` activo (sin rutas
fijas) e instancia el entorno automáticamente la primera vez.

Servidor Pluto general (eliges el notebook desde el navegador):

```bash
julia --project=. scripts/pluto.jl        # o: make pluto
```

Abrir directamente un notebook concreto:

```bash
julia --project=. scripts/pluto.jl bgw          # BGW básico             (make bgw)
julia --project=. scripts/pluto.jl parejas      # Parejas básico         (make parejas)
julia --project=. scripts/pluto.jl bgw-dane     # BGW con datos DANE     (make bgw-dane)
julia --project=. scripts/pluto.jl parejas-dane # Árbol genealógico DANE (make parejas-dane)
```

Los notebooks se organizan en dos carpetas: `notebooks/basic/` ilustra cómo
funciona cada modelo de forma interactiva, y `notebooks/results/` aplica los
modelos a datos reales del Censo DANE 2018. Todos cargan el código únicamente de
`src/` (no redefinen los modelos). El modelo por parejas está documentado en
[`docs/modelo-parejas.md`](docs/modelo-parejas.md).

## Catálogo reproducible de experimentos (sin Pluto)

Para resultados citables sin abrir Pluto, `scripts/experimentos_parejas.jl` corre
un catálogo fijo de experimentos (completo, Romeo y Julieta, cadena con puente,
estrella, anillo, grid, bipartito, tripartito y componentes) con semillas fijas y
escribe a disco en `output/`:

- Un **PNG** del grafo inicial, del conteo por apellido, de la diversidad y de la
  distribución de supervivientes, por cada experimento.
- Un **CSV** (`output/resumen_experimentos.csv`) con el resumen agregado de las
  repeticiones: probabilidad de absorción, nº medio de supervivientes,
  probabilidad de ganar por apellido, probabilidad de que el puente muera primero
  y mutaciones medias.

```bash
julia --project=. scripts/experimentos_parejas.jl   # o: make experimentos
```

Usa solo `DelimitedFiles` (librería estándar), sin dependencias nuevas. Se
parametriza por variables de entorno (con valores por defecto reproducibles):

| Variable | Qué controla | Default |
| --- | --- | --- |
| `EXP_N_REPS` | repeticiones por experimento | `200` |
| `EXP_PAREJAS` | parejas por apellido en el estado inicial | `3` |
| `EXP_SEMILLA_BASE` | semilla base de las repeticiones | `1` |
| `EXP_MAX_GEN` | tope de generaciones por corrida | `1000` |
| `EXP_SEMILLA_FIG` | semilla de la corrida de ejemplo (figuras) | `7` |
| `EXP_OUTDIR` | carpeta de salida | `output` |

Ejemplo con otra configuración:

```bash
EXP_N_REPS=500 EXP_PAREJAS=5 julia --project=. scripts/experimentos_parejas.jl
```

La carpeta `output/` está en `.gitignore`: es salida generada, se regenera con el
script.

## Inicializar el proyecto

El proyecto usa su propio entorno (`Project.toml`), independiente de tu entorno
global de Julia. Esto evita conflictos de dependencias. **Solo hay que hacerlo
una vez.**

Desde la raíz del proyecto, instala y precompila las dependencias:

```powershell
julia --project=. -e "import Pkg; Pkg.instantiate()"
```

Si aparece un error de resolución de paquetes (por ejemplo con `Qt6Base_jll`),
actualiza el registro y vuelve a intentarlo:

```powershell
julia --project=. -e "import Pkg; Pkg.Registry.update(); Pkg.resolve(); Pkg.instantiate()"
```

> **Importante:** usa siempre la bandera `--project=.` al ejecutar cualquier
> comando. Sin ella, Julia usaría tu entorno global en lugar del entorno
> aislado del proyecto.

## Ejecutar los tests

Todos los comandos se ejecutan **desde la raíz del proyecto**.

Tests de los métodos numéricos (bisección, Newton, secante):

```powershell
julia --project=. test/runtests.jl
```

Tests del modelo de Markov y la simulación de extinción:

```powershell
julia --project=. test/test_markov.jl
```

El segundo imprime el progreso en vivo (bloques `1/7`, `2/7`, ...) con marcas de
tiempo, para que puedas ver que está avanzando. La primera ejecución tarda más
porque Julia precompila las dependencias; las siguientes son más rápidas.

Salida esperada (resumen):

```
Metodos numericos |   25     25   (todos pasan)
Modelo de Markov / Galton-Watson |   28     28   (todos pasan)
```

## Ejecutar la simulación

```powershell
julia --project=. scripts/run_simulation.jl
```

Muestra, para un modelo de natalidad alta y otro de natalidad baja:

- La media de descendencia `m` y el régimen (subcrítico / crítico / supercrítico).
- La probabilidad de extinción teórica `q`.
- Trayectorias de ejemplo de linajes individuales.
- Una validación Monte Carlo que compara la frecuencia de extinción simulada
  contra la `q` teórica.

## Benchmark

Para medir el rendimiento de los métodos numéricos:

```powershell
julia --project=. test/benchmarks.jl
```

Usa el paquete `BenchmarkTools` si está disponible; si no, recurre a una
medición simple con la librería estándar. Para instalar `BenchmarkTools` en el
entorno del proyecto:

```powershell
julia --project=. -e "import Pkg; Pkg.add(\"BenchmarkTools\")"
```

## Solución de problemas

- **`The active manifest file has dependencies that were resolved with a
  different julia version` seguido de fallos de precompilación de `Plots` /
  `GR` / `GR_jll` / `Libtiff_jll` (p. ej. `Zstd_jll ... is required but does
  not seem to be installed`)**: tu versión de Julia no coincide con la del
  `Manifest.toml`. Asegúrate de usar Julia 1.12 en esta carpeta (ver *Versión
  de Julia fija en todos los dispositivos*):

  ```bash
  juliaup override set 1.12
  julia --version            # debe decir 1.12.x
  ```

  Si el `Manifest.toml` quedó resuelto con otra versión, regénéralo con la
  versión correcta:

  ```bash
  rm Manifest.toml
  julia --project=. -e "import Pkg; Pkg.resolve(); Pkg.instantiate()"
  ```

- **`Package X not found in current path`**: olvidaste `--project=.`, o falta
  ejecutar `Pkg.instantiate()`. Revisa la sección *Inicializar el proyecto*.
- **El test del modelo de Markov parece colgado**: no lo está; el bloque de
  Monte Carlo tarda unos segundos. Verás el mensaje `7/7 ... puede tardar unos
  seg` justo antes. La primera ejecución además precompila `Distributions`.
- **Error de resolución con `Qt6Base_jll` u otro `_jll`**: proviene del entorno
  global. El `Project.toml` del proyecto evita este problema; asegúrate de usar
  `--project=.`.
