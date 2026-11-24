# src/core/visualization.jl
#
# Utilidades de graficacion para el modelo por parejas: dibuja el grafo de
# apellidos y grafica conteos, diversidad y estadisticas de supervivencia a
# partir de las salidas del modelo. No contiene logica de modelo.
#
# Documentacion completa (graficas disponibles, uso): docs/code/visualization.md
#
# Usa Plots.jl y, para dibujar los grafos, Graphs.jl + GraphRecipes.jl.
# Requiere que experimentMetrics.jl este cargado antes que este archivo.

using Plots
using Graphs
using GraphRecipes
using StatsBase

"""
    _mapa_indices(compatibles) -> (ids, idx)

Devuelve el vector ordenado de ids de apellido `ids` y un Dict `idx` que mapea
cada id a su posicion 1..n (indice de vertice para Graphs.jl).
"""
function _mapa_indices(compatibles::Dict{Int,Vector{Int}})
    ids = sort(collect(keys(compatibles)))
    idx = Dict(a => k for (k, a) in enumerate(ids))
    return ids, idx
end

"""
    dibujar_grafo(compatibles; conteos=nothing, titulo="Grafo de apellidos",
                  mostrar_autolazo=false)

Dibuja el grafo dirigido de apellidos. Cada nodo es un apellido; su etiqueta es
el id (y el conteo si se pasa `conteos`). El tamano del nodo es proporcional al
conteo cuando se proporciona.

- `compatibles`: la adyacencia de la topologia.
- `conteos`: Dict apellido => personas (opcional, para tamano/etiqueta).
- `mostrar_autolazo`: si `false` (por defecto) no dibuja los auto-lazos i->i
  para no recargar la figura.
"""
function dibujar_grafo(compatibles::Dict{Int,Vector{Int}};
                       conteos::Union{Dict{Int,Int},Nothing}=nothing,
                       titulo::AbstractString="Grafo de apellidos",
                       mostrar_autolazo::Bool=false)
    ids, idx = _mapa_indices(compatibles)
    n = length(ids)
    g = SimpleDiGraph(n)
    n_aristas = 0
    for i in ids
        for j in compatibles[i]
            if (i != j || mostrar_autolazo) && haskey(idx, j)
                add_edge!(g, idx[i], idx[j])
                n_aristas += 1
            end
        end
    end

    # Caso borde: graphplot no puede inferir layout de un grafo con <=1 nodo o
    # sin aristas (p.ej. tras la fijacion a un solo apellido). Dibujamos un
    # marcador de posicion claro en vez de fallar.
    if n <= 1 || n_aristas == 0
        etiqueta = if n == 0
            "Sin apellidos"
        else
            txt = join([conteos === nothing ? string(a) : string(a, " (", get(conteos, a, 0), ")") for a in ids], ", ")
            "Apellido(s) sin aristas: " * txt
        end
        plt = plot(; title = titulo, legend = false, framestyle = :none,
                   xlims = (0, 1), ylims = (0, 1))
        annotate!(plt, 0.5, 0.5, text(etiqueta, 10))
        return plt
    end

    etiquetas = if conteos === nothing
        [string(a) for a in ids]
    else
        [string(a, " (", get(conteos, a, 0), ")") for a in ids]
    end

    return graphplot(g;
                     names = etiquetas,
                     nodeshape = :circle,
                     nodesize = 0.15,
                     nodecolor = :lightblue,
                     curves = false,
                     arrow = true,
                     title = titulo,
                     fontsize = 8)
end

"""
    graficar_conteos(historial; apellidos=nothing, titulo="Conteo por apellido")

Grafica el conteo de cada apellido a lo largo de las generaciones (una curva por
apellido). Si no se pasan `apellidos`, usa todos los presentes en la generacion
inicial.
"""
function graficar_conteos(historial::Vector{Vector{Int}};
                          apellidos::Union{Vector{Int},Nothing}=nothing,
                          titulo::AbstractString="Conteo por apellido")
    ap = apellidos === nothing ? sort(unique(historial[1])) : apellidos
    M = conteos_por_gen(historial, ap)
    gens = 0:(size(M, 1) - 1)
    plt = plot(; title = titulo, xlabel = "Generacion", ylabel = "Personas",
               legend = :outertopright)
    for (k, a) in enumerate(ap)
        plot!(plt, gens, M[:, k]; label = "Apellido $a", lw = 2)
    end
    return plt
end

"""
    graficar_diversidad(historial; titulo="Apellidos vivos por generacion")

Grafica el numero de apellidos distintos vivos en cada generacion (serie no
creciente que termina en 1 si hubo fijacion).
"""
function graficar_diversidad(historial::Vector{Vector{Int}};
                             titulo::AbstractString="Apellidos vivos por generacion")
    vivos = apellidos_vivos_por_gen(historial)
    gens = 0:(length(vivos) - 1)
    return plot(gens, vivos;
                title = titulo, xlabel = "Generacion",
                ylabel = "Apellidos vivos", lw = 2, legend = false,
                marker = :circle, markersize = 2)
end

"""
    graficar_dist_n_vivos(rep; titulo="Distribucion de apellidos supervivientes")

Grafico de barras con la fraccion de corridas que terminaron con k apellidos
vivos, tomado de la salida de `correr_repeticiones`. Muestra si el sistema
tiende a fijar uno solo o a mantener coexistencia.
"""
function graficar_dist_n_vivos(rep; titulo::AbstractString="Distribucion de apellidos supervivientes")
    ks = sort(collect(keys(rep.dist_n_vivos)))
    fracs = [rep.dist_n_vivos[k] for k in ks]
    return bar(string.(ks), fracs;
               title = titulo, xlabel = "Apellidos vivos al final",
               ylabel = "Fraccion de corridas", legend = false, ylims = (0, 1))
end

"""
    graficar_prob_ganar(rep; titulo="Probabilidad de ser unico superviviente")

Grafico de barras con la probabilidad de que cada apellido sea el unico
superviviente (solo corridas que fijaron un apellido).
"""
function graficar_prob_ganar(rep; titulo::AbstractString="Probabilidad de ser unico superviviente")
    ap = sort(collect(keys(rep.prob_ganar)))
    probs = [rep.prob_ganar[a] for a in ap]
    return bar(string.(ap), probs;
               title = titulo, xlabel = "Apellido", ylabel = "P(unico superviviente)",
               legend = false, ylims = (0, 1))
end
