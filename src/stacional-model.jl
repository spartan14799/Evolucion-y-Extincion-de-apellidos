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
Simula la evolución hasta la fijación de un único apellido, o hasta alcanzar el límite de generaciones.
"""
function simular_fijacion(estado_inicial::Vector{Int}; max_generaciones::Int=10_000, rng::AbstractRNG=Random.default_rng())
    poblacion = copy(estado_inicial)
    historial = [copy(poblacion)]
    generaciones = 0
    
    # Continúa mientras haya más de 1 apellido Y no choquemos con el límite
    while length(unique(poblacion)) > 1 && generaciones < max_generaciones
        poblacion = siguiente_generacion(poblacion, rng)
        push!(historial, copy(poblacion))
        generaciones += 1
    end
    
    # Determinamos si el ciclo terminó porque ganó alguien, o porque se nos acabó el tiempo
    logro_fijacion = length(unique(poblacion)) == 1
    ganador = logro_fijacion ? poblacion[1] : nothing
    
    return (ganador = ganador, generaciones = generaciones, historial = historial, fijacion = logro_fijacion)
end
