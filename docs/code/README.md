# Documentación del código (`docs/code/`)

Una ficha por archivo de código del proyecto: qué hace, cómo lo hace, qué
funciones/tipos expone (con firma), de qué depende y cómo se usa. El encabezado
de cada archivo `.jl` sólo lleva una descripción breve y apunta a la ficha
correspondiente de esta carpeta.

## Arquitectura en una vista

```
src/
├── bgwModel.jl          MODELO 1: Bienaymé-Galton-Watson (extinción de linajes)
├── coupleModel.jl       MODELO 2: modelo por parejas sobre el grafo de apellidos
├── core/                requisitos genéricos de los modelos
│   ├── numericalMethods.jl   búsqueda de raíces (bisección, Newton, secante)
│   ├── graphTopologies.jl    generadores de grafos de compatibilidad
│   ├── experimentMetrics.jl  métricas y runner del modelo por parejas
│   └── visualization.jl      gráficas (grafo, conteos, diversidad, barras)
└── experiments/         datos reales y orquestación para los notebooks de results/
    ├── daneData.jl           datos del Censo DANE 2018 (solo datos)
    ├── daneBgwRunners.jl     Monte Carlo BGW sobre los datos DANE
    └── daneCoupleRunners.jl  árbol genealógico por parejas sobre los datos DANE
```

Regla de dependencias: `experiments/` depende de los modelos y de `core/`;
`core/` da soporte a los modelos; los modelos no dependen de `core/` salvo BGW,
que usa `core/numericalMethods.jl`. Nada en `src/` depende de los notebooks.

## Fichas

### Modelos (`src/`)
- [bgwModel.md](bgwModel.md) — modelo de ramificación Bienaymé-Galton-Watson.
- [coupleModel.md](coupleModel.md) — modelo por parejas sobre el grafo.

### Soporte (`src/core/`)
- [numericalMethods.md](numericalMethods.md) — bisección, Newton, secante.
- [graphTopologies.md](graphTopologies.md) — generadores de topologías del grafo.
- [experimentMetrics.md](experimentMetrics.md) — métricas y runner de parejas.
- [visualization.md](visualization.md) — utilidades de graficación.

### Experimentos con datos reales (`src/experiments/`)
- [daneData.md](daneData.md) — datos del Censo DANE 2018.
- [daneBgwRunners.md](daneBgwRunners.md) — orquestación Monte Carlo BGW.
- [daneCoupleRunners.md](daneCoupleRunners.md) — orquestación del árbol por parejas.

### Scripts (`scripts/`)
- [scripts.md](scripts.md) — `pluto.jl`, `run_simulation.jl`, `experimentos_parejas.jl`.

### Notebooks (`notebooks/`)
- [notebooks.md](notebooks.md) — `basic/` (ilustración) y `results/` (datos reales).

### Tests (`test/`)
- [tests.md](tests.md) — `runtests.jl`, `test_markov.jl`, `benchmarks.jl`.
