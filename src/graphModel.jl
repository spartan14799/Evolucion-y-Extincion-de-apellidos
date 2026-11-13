# src/graphModel.jl

struct GrafoApellidos
    apellidos::Vector{Int}
    conteos::Dict{Int,Int}
    compatibles::Dict{Int,Vector{Int}}
    probabilidades::Dict{Int,Dict{Int,Float64}}
end

"""
    contar_apellidos(poblacion::Vector{Int}) -> Dict{Int,Int}

Cuenta la cantidad de personas por cada ID de apellido en la generación actual.
"""
function contar_apellidos(poblacion::Vector{Int})
    conteos = Dict{Int,Int}()
    for id in poblacion
        conteos[id] = get(conteos, id, 0) + 1
    end
    return conteos
end

"""
    GrafoApellidos(poblacion; compatibles=nothing)

Constructor que inicializa los vértices y conteos. Si no se provee `compatibles`,
se asume un grafo completo (cada apellido es compatible con todos, incluido sí mismo).
"""
function GrafoApellidos(poblacion::Vector{Int}; compatibles=nothing)
    conteos = contar_apellidos(poblacion)
    apellidos = sort(collect(keys(conteos)))

    comp = Dict{Int,Vector{Int}}()
    if compatibles === nothing
        # Por defecto: grafo completo
        for i in apellidos
            comp[i] = copy(apellidos)
        end
    else
        comp = compatibles
    end

    # Las probabilidades inician vacías; se calcularán en la Task 2.
    prob = Dict{Int,Dict{Int,Float64}}()

    return GrafoApellidos(apellidos, conteos, comp, prob)
end

"""
    probabilidades_proporcionales(grafo::GrafoApellidos) -> Dict{Int,Dict{Int,Float64}}

Calcula las probabilidades de transición P(i → j) de forma proporcional
al conteo del apellido destino j, considerando solo los apellidos
compatibles con i.
"""
function probabilidades_proporcionales(grafo::GrafoApellidos)
    probabilidades = Dict{Int,Dict{Int,Float64}}()
    
    for i in grafo.apellidos
        probabilidades[i] = Dict{Int,Float64}()
        # Denominador: suma de personas de los apellidos compatibles con i
        total_compatibles = 0
        if haskey(grafo.compatibles, i)
            for j in grafo.compatibles[i]
                total_compatibles += get(grafo.conteos, j, 0)
            end
        end
        
        if total_compatibles > 0
            for j in grafo.compatibles[i]
                conteo_j = get(grafo.conteos, j, 0)
                probabilidades[i][j] = conteo_j / total_compatibles
            end
        end
    end
    
    return probabilidades
end

"""
    probabilidades_salientes(grafo::GrafoApellidos, i::Int) -> Dict{Int,Float64}

Devuelve el diccionario de probabilidades salientes para el vértice `i`.
"""
function probabilidades_salientes(grafo::GrafoApellidos, i::Int)
    return get(grafo.probabilidades, i, Dict{Int,Float64}())
end
