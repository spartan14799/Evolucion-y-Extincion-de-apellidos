using Distributions

"""
    ModeloPBGW

Estructura que representa un Proceso de Bienayme-Galton-Watson para unico apellidos.
Almacena la distribucion de probabilidad de descendencia.
"""
struct ModeloPBGW
    pk::Vector{Float64}

    function ModeloPBGW(pk::Vector{Float64})
        if any(p -> p < 0.0, pk) error("Probabilidades negativas.") end
        if !isapprox(sum(pk), 1.0, atol=1e-5) error("No suman 1.0") end
        new(pk)
    end

    function ModeloPBGW(tasa_natalidad::Float64; max_hijos::Int=15)
        # Poisson lambda = tasa_natalidad
        dist = Poisson(tasa_natalidad)
        
        pk = [pdf(dist, k) for k in 0:max_hijos]
        
        # Normalizamos para segurar que sumen exactamente 1
        pk = pk ./ sum(pk)
    
        new(pk)
    end
end

"""
Calcula el promedio de hijos por individuo (Esperanza m = E[X]).
"""
function esperanza(modelo::ModeloPBGW)
    return sum((k - 1) * modelo.pk[k] for k in 1:length(modelo.pk))
end

"""
Verifica si la extincion es matematicamente 100% segura (m <= 1).
"""
function es_extincion_segura(modelo::ModeloPBGW)
    return esperanza(modelo) <= 1.0
end

"""
Evalua la Funcion Generadora de Probabilidad
"""
function fgp(modelo::ModeloPBGW, t::Float64)
    return sum(modelo.pk[k] * (t^(k - 1)) for k in 1:length(modelo.pk))
end

