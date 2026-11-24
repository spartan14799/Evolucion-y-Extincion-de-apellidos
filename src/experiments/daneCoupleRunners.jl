# src/experiments/daneCoupleRunners.jl
#
# Orquestacion del MODELO POR PAREJAS aplicada a los datos del Censo DANE 2018
# (arbol genealogico). Alimenta el notebook notebooks/results/couple_genealogy_dane.jl:
# corre la simulacion con el motor `siguiente_generacion_parejas` de
# src/coupleModel.jl y devuelve el historial y los metadatos para que el notebook
# solo dibuje. No reimplementa el modelo.
#
# Documentacion completa (funciones, uso): docs/code/daneCoupleRunners.md
#
# Requiere cargados antes: src/coupleModel.jl, src/core/graphTopologies.jl y
# src/experiments/daneData.jl.

using Random
using StatsBase

"""
    simular_arbol_dane(; modo_inicial=:equitativo, parejas=10,
                       modo_prob=:ponderado, semilla=7, max_generaciones=25)
        -> NamedTuple

Simula el arbol genealogico por parejas sobre los 10 apellidos del Censo DANE
2018. Devuelve la informacion necesaria para que el notebook solo dibuje:
- `historial`   : Vector con la poblacion de cada generacion (indice 1 = gen 0).
- `generaciones`: numero de transiciones simuladas.
- `caidos`      : Dict apellido_id => generacion en que se extinguio.
- `ganador`     : id del apellido dominante al final.
- `fijado`      : true si al final queda un unico apellido.
- `n_total`     : tamano de la poblacion (2 * parejas).

Modos de poblacion inicial:
- `:equitativo`: reparte las parejas ciclicamente entre los 10 apellidos.
- `:dane`      : muestrea el apellido de cada pareja proporcional a `pesos_dane()`.
"""
function simular_arbol_dane(; modo_inicial::Symbol = :equitativo,
                            parejas::Int = 10,
                            modo_prob::Symbol = :ponderado,
                            semilla::Int = 7,
                            max_generaciones::Int = 25)
    n_total = 2 * parejas
    rng = MersenneTwister(semilla)

    apellidos_pareja = if modo_inicial == :equitativo
        [mod1(k, 10) for k in 1:parejas]
    else
        sample(rng, 1:10, Weights(pesos_dane()), parejas)
    end

    poblacion = Int[]
    for a in apellidos_pareja
        append!(poblacion, (a, a))
    end

    # Grafo completo sobre los 10 apellidos (todos pueden emparejarse con todos).
    compatibles = topo_completo(10)

    historial = [copy(poblacion)]
    caidos = Dict{Int,Int}()
    for g in 1:max_generaciones
        actual = historial[g]
        hijos, _ = siguiente_generacion_parejas(actual, compatibles, rng;
                                                modo_prob = modo_prob)
        push!(historial, hijos)
        vivos = Set(hijos)
        for a in unique(historial[1])
            if a ∉ vivos && !haskey(caidos, a)
                caidos[a] = g
            end
        end
        length(vivos) == 1 && break
    end

    generaciones = length(historial) - 1
    ganador = historial[end][1]
    fijado = length(unique(historial[end])) == 1

    return (historial = historial,
            generaciones = generaciones,
            caidos = caidos,
            ganador = ganador,
            fijado = fijado,
            n_total = n_total)
end
