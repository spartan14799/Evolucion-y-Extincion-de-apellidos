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