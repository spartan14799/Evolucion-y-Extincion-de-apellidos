using Random

"""
Simula una transición de estado en la población bajo el modelo de parejas de tamaño constante.
"""
function siguiente_generacion(poblacion::Vector{Int}, rng::AbstractRNG=Random.default_rng())
    mezcla = shuffle(rng, poblacion)
    nueva_poblacion = Int[]
    
    for i in 1:2:length(mezcla)
        apellido_ganador = mezcla[i]
        push!(nueva_poblacion, apellido_ganador, apellido_ganador)
    end
    
    return nueva_poblacion
end

"""
Simula la evolución hasta la fijación de un único apellido, partiendo de un estado inicial dado.
"""
function simular_fijacion(estado_inicial::Vector{Int}; rng::AbstractRNG=Random.default_rng())
    poblacion = copy(estado_inicial)
    historial = [copy(poblacion)]
    generaciones = 0
    
    while length(unique(poblacion)) > 1
        poblacion = siguiente_generacion(poblacion, rng)
        push!(historial, copy(poblacion))
        generaciones += 1
    end
    
    return (ganador = poblacion[1], generaciones = generaciones, historial = historial)
end
