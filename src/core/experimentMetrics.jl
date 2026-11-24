# src/core/experimentMetrics.jl
#
# Metricas y runner de experimentos sobre el modelo por parejas
# (src/coupleModel.jl). Post-procesa el `historial` y las parejas forzadas que
# devuelve `simular_parejas` (diversidad, orden de extincion, mutaciones) y
# orquesta corridas repetidas con estadisticas agregadas. No toca el motor.
#
# Documentacion completa (metricas, runner, uso): docs/code/experimentMetrics.md
#
# Requiere que coupleModel.jl este cargado antes que este archivo.

using Random
using StatsBase

"""
    apellidos_vivos_por_gen(historial) -> Vector{Int}

Numero de apellidos distintos vivos en cada generacion. Serie no creciente.
"""
function apellidos_vivos_por_gen(historial::Vector{Vector{Int}})
    return [length(unique(gen)) for gen in historial]
end

"""
    conteos_por_gen(historial, apellidos) -> Matrix{Int}

Matriz (generaciones+1) x length(apellidos) con el conteo de cada apellido en
cada generacion. `apellidos` fija el orden de columnas.
"""
function conteos_por_gen(historial::Vector{Vector{Int}}, apellidos::Vector{Int})
    G = length(historial)
    M = zeros(Int, G, length(apellidos))
    col = Dict(a => k for (k, a) in enumerate(apellidos))
    for (g, gen) in enumerate(historial)
        for (a, c) in countmap(gen)
            haskey(col, a) && (M[g, col[a]] = c)
        end
    end
    return M
end

"""
    generacion_de_extincion(historial) -> Dict{Int,Int}

Para cada apellido presente en la generacion 0, la generacion (base 0) en que su
conteo llega a cero por primera vez. Si nunca se extingue, `typemax(Int)`.
"""
function generacion_de_extincion(historial::Vector{Vector{Int}})
    iniciales = unique(historial[1])
    extincion = Dict{Int,Int}(a => typemax(Int) for a in iniciales)
    for (g, gen) in enumerate(historial)
        vivos = Set(gen)
        for a in iniciales
            if extincion[a] == typemax(Int) && !(a in vivos)
                extincion[a] = g - 1
            end
        end
    end
    return extincion
end

"""
    orden_de_extincion(historial) -> Vector{Pair{Int,Int}}

Apellidos ordenados por la generacion de su extincion (primero el que muere
antes). `apellido => generacion`; los supervivientes quedan al final con
`typemax(Int)`.
"""
function orden_de_extincion(historial::Vector{Vector{Int}})
    ext = generacion_de_extincion(historial)
    return sort(collect(ext); by = p -> p.second)
end

"""
    total_mutaciones(resultado::ResultadoParejas) -> Int

Numero total de emparejamientos forzados (mutaciones definidas) a lo largo de
toda la simulacion.
"""
total_mutaciones(resultado::ResultadoParejas) = sum(length, resultado.forzadas_por_gen; init = 0)

"""
    mutaciones_por_gen(resultado::ResultadoParejas) -> Vector{Int}

Numero de emparejamientos forzados en cada generacion.
"""
mutaciones_por_gen(resultado::ResultadoParejas) = [length(f) for f in resultado.forzadas_por_gen]

"""
    resumen_experimento(resultado::ResultadoParejas; puente=nothing) -> NamedTuple

Resume una corrida de `simular_parejas`.

Campos:
- `generaciones`, `absorbido`, `transitorio`: del resultado.
- `apellidos_iniciales`, `vivos_final`.
- `n_vivos_final`: cuantos apellidos sobreviven.
- `orden_extincion`.
- `mutaciones`: total de emparejamientos forzados.
- `puente_muere_primero`: `true/false/nothing`. `nothing` si no se paso `puente`
  o si el puente no se extingue. `true` si el apellido `puente` fue el primero en
  extinguirse entre los que se extinguen.
"""
function resumen_experimento(resultado::ResultadoParejas; puente::Union{Int,Nothing}=nothing)
    historial = resultado.historial
    orden = orden_de_extincion(historial)
    iniciales = sort(unique(historial[1]))
    vivos_final = sort(unique(historial[end]))

    puente_muere_primero = nothing
    if puente !== nothing
        extintos = [p.first for p in orden if p.second != typemax(Int)]
        puente_muere_primero = isempty(extintos) ? false : (first(extintos) == puente)
    end

    return (generaciones = resultado.generaciones,
            absorbido = resultado.absorbido,
            transitorio = resultado.transitorio,
            apellidos_iniciales = iniciales,
            vivos_final = vivos_final,
            n_vivos_final = length(vivos_final),
            orden_extincion = orden,
            mutaciones = total_mutaciones(resultado),
            puente_muere_primero = puente_muere_primero)
end

"""
    poblacion_desde_parejas(parejas_por_apellido::Dict{Int,Int}) -> Vector{Int}

Construye el estado inicial a partir del NUMERO DE PAREJAS de cada apellido.
Cada pareja son 2 personas, de modo que la poblacion total es siempre par.
Ej: `Dict(1 => 2, 2 => 3)` -> apellido 1 con 4 personas, apellido 2 con 6.
"""
function poblacion_desde_parejas(parejas_por_apellido::Dict{Int,Int})
    estado = Int[]
    for a in sort(collect(keys(parejas_por_apellido)))
        append!(estado, fill(a, 2 * parejas_por_apellido[a]))
    end
    return estado
end

"""
    correr_experimento(estado_inicial, compatibles; modo_prob=:uniforme_aristas,
                       max_generaciones=1000, semilla=nothing) -> ResultadoParejas

Atajo para correr una sola simulacion del modelo por parejas.
"""
function correr_experimento(estado_inicial::Vector{Int}, compatibles::Dict{Int,Vector{Int}};
                            modo_prob::Symbol = :uniforme_aristas,
                            max_generaciones::Int = 1000,
                            semilla::Union{Int,Nothing} = nothing)
    parametros = ParametrosParejas(max_generaciones = max_generaciones, modo_prob = modo_prob)
    rng = semilla === nothing ? Random.default_rng() : MersenneTwister(semilla)
    return simular_parejas(estado_inicial, compatibles; parametros = parametros, rng = rng)
end

"""
    correr_repeticiones(estado_inicial, compatibles; modo_prob=:uniforme_aristas,
                        max_generaciones=1000, n_reps=100, semilla_base=1,
                        puente=nothing) -> NamedTuple

Corre el experimento `n_reps` veces con semillas reproducibles
(`semilla_base .+ (0:n_reps-1)`) y agrega estadisticas:

- `n_reps`.
- `prob_absorcion`: fraccion de corridas que alcanzaron un estado absorbente
  (fijacion o equilibrio permanente) dentro del tope.
- `prob_transitorio`: fraccion que termino en estado transitorio.
- `gen_media`, `gen_std`: generaciones hasta parar.
- `dist_n_vivos`: Dict k => fraccion de corridas que terminaron con k apellidos.
- `prob_ganar`: Dict apellido => fraccion de corridas en que fue el UNICO
  superviviente (solo cuenta cuando termina con 1 apellido).
- `mutaciones_media`: promedio de emparejamientos forzados por corrida.
- `prob_puente_primero`: fraccion de corridas en que `puente` se extingue primero.
- `resumenes`: vector con el resumen de cada corrida.
"""
function correr_repeticiones(estado_inicial::Vector{Int}, compatibles::Dict{Int,Vector{Int}};
                             modo_prob::Symbol = :uniforme_aristas,
                             max_generaciones::Int = 1000,
                             n_reps::Int = 100,
                             semilla_base::Int = 1,
                             puente::Union{Int,Nothing} = nothing)
    resumenes = Vector{Any}(undef, n_reps)
    generaciones = zeros(Int, n_reps)
    absorciones = 0
    transitorios = 0
    ganadores = Int[]
    n_vivos = Int[]
    mutaciones = zeros(Int, n_reps)
    puente_primero = 0
    puente_contables = 0

    for r in 1:n_reps
        res = correr_experimento(estado_inicial, compatibles;
                                 modo_prob = modo_prob,
                                 max_generaciones = max_generaciones,
                                 semilla = semilla_base + r - 1)
        resumen = resumen_experimento(res; puente = puente)
        resumenes[r] = resumen
        generaciones[r] = res.generaciones
        mutaciones[r] = resumen.mutaciones
        push!(n_vivos, resumen.n_vivos_final)
        res.absorbido && (absorciones += 1)
        res.transitorio && (transitorios += 1)
        if resumen.n_vivos_final == 1
            push!(ganadores, res.apellidos_final[1])
        end
        if resumen.puente_muere_primero !== nothing
            puente_contables += 1
            resumen.puente_muere_primero && (puente_primero += 1)
        end
    end

    prob_ganar = Dict{Int,Float64}()
    for (a, c) in countmap(ganadores)
        prob_ganar[a] = c / n_reps
    end

    dist_n_vivos = Dict{Int,Float64}()
    for (k, c) in countmap(n_vivos)
        dist_n_vivos[k] = c / n_reps
    end

    return (n_reps = n_reps,
            prob_absorcion = absorciones / n_reps,
            prob_transitorio = transitorios / n_reps,
            gen_media = mean(generaciones),
            gen_std = n_reps > 1 ? std(generaciones) : 0.0,
            dist_n_vivos = dist_n_vivos,
            prob_ganar = prob_ganar,
            mutaciones_media = mean(mutaciones),
            prob_puente_primero = puente_contables > 0 ? puente_primero / puente_contables : nothing,
            resumenes = resumenes)
end
