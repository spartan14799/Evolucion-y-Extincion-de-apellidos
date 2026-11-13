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
