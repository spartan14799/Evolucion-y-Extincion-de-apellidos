# src/stacional-model.jl

using Random

# Incluir el motor de grafos
include("graphModel.jl")

"""
    ParametrosModelo

Parámetros para la simulación del modelo estacionario.
`fijar_probabilidades` permite inyectar una función que calcule probabilidades custom 
en cada generación. Si es `nothing`, se usa el default proporcional del grafo.
"""
Base.@kwdef struct ParametrosModelo
    max_generaciones::Int = 1000
    fijar_probabilidades::Union{Function, Nothing} = nothing
end

"""
    simular_fijacion(estado_inicial::Vector{Int}; parametros=ParametrosModelo(), rng=Random.default_rng())

Simula la evolución de la población hasta que todos los individuos compartan el mismo apellido 
(fijación) o se alcance el límite máximo de generaciones.

Retorna una NamedTuple: `(ganador, generaciones, historial, fijacion)`
"""
function simular_fijacion(estado_inicial::Vector{Int}; 
                          parametros::ParametrosModelo = ParametrosModelo(), 
                          rng::AbstractRNG = Random.default_rng())
    
    historial = Vector{Vector{Int}}()
    push!(historial, copy(estado_inicial))
    
    poblacion_actual = copy(estado_inicial)
    generaciones = 0
    fijacion = false
    ganador = nothing
    
    # Chequeo de estado absorbente inicial
    if length(unique(poblacion_actual)) <= 1
        fijacion = true
        ganador = isempty(poblacion_actual) ? nothing : poblacion_actual[1]
        return (ganador=ganador, generaciones=generaciones, historial=historial, fijacion=fijacion)
    end
    
    for g in 1:parametros.max_generaciones
        # Si el usuario pasó una estrategia custom de cálculo, la ejecutamos
        probs_custom = parametros.fijar_probabilidades !== nothing ? parametros.fijar_probabilidades(poblacion_actual) : nothing
        
        # Transición usando el motor estocástico basado en el grafo dirigido
        poblacion_actual = siguiente_generacion_grafo(poblacion_actual, rng; probabilidades_custom=probs_custom)
        
        push!(historial, copy(poblacion_actual))
        generaciones += 1
        
        # Verificar estado absorbente (fijación)
        unicos = unique(poblacion_actual)
        if length(unicos) == 1
            fijacion = true
            ganador = unicos[1]
            break
        end
    end
    
    return (ganador=ganador, generaciones=generaciones, historial=historial, fijacion=fijacion)
end
