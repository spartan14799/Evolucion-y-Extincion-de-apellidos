# src/coupleModel.jl
#
# MODELO 2: modelo por parejas sobre el grafo de apellidos (modelo vigente).
# Nodo = apellido; las aristas (`compatibles`) definen que apellidos pueden
# emparejarse. La unidad reproductiva es la PAREJA y cada una aporta 2 hijos,
# de modo que N = 2P se conserva. Motor independiente: simula la evolucion
# hasta fijacion, equilibrio permanente o tope de generaciones.
#
# Documentacion completa (reglas, funciones, uso): docs/code/coupleModel.md

using Random
using StatsBase

"""
    contar(poblacion::Vector{Int}) -> Dict{Int,Int}

Conteo de personas por apellido en la poblacion actual.
"""
function contar(poblacion::Vector{Int})
    d = Dict{Int,Int}()
    for a in poblacion
        d[a] = get(d, a, 0) + 1
    end
    return d
end

"""
    ParametrosParejas

Parametros del modelo por parejas.

- `max_generaciones`: tope de generaciones si no se alcanza un estado absorbente.
- `modo_prob`: `:uniforme_aristas` (default) o `:ponderado`.
    - `:uniforme_aristas`: cada persona de apellido `i` elige con quien
      emparejarse de forma uniforme entre las aristas salientes de `i` cuyos
      apellidos siguen VIVOS. Si un apellido se extingue, su arista desaparece y
      la probabilidad se reparte entre las aristas restantes.
    - `:ponderado`: la probabilidad de emparejarse con `j` es proporcional al
      conteo de `j` (entre los compatibles de `i`), calculado una vez por
      generacion.
"""
Base.@kwdef struct ParametrosParejas
    max_generaciones::Int = 1000
    modo_prob::Symbol = :uniforme_aristas
end

"""
    pesos_emparejamiento(apellido_i, compatibles, conteos, vivos, modo) -> (destinos, pesos)

Devuelve los apellidos destino posibles para una persona de apellido `apellido_i`
y sus pesos relativos, segun el modo. Solo se consideran apellidos VIVOS
(presentes en la generacion). Los pesos no tienen que estar normalizados.
"""
function pesos_emparejamiento(apellido_i::Int,
                              compatibles::Dict{Int,Vector{Int}},
                              conteos::Dict{Int,Int},
                              vivos::Set{Int},
                              modo::Symbol)
    vecinos = haskey(compatibles, apellido_i) ? compatibles[apellido_i] : Int[]
    destinos = [j for j in vecinos if j in vivos]
    if isempty(destinos)
        return (Int[], Float64[])
    end
    if modo == :uniforme_aristas
        pesos = fill(1.0, length(destinos))
    elseif modo == :ponderado
        pesos = [Float64(get(conteos, j, 0)) for j in destinos]
        if sum(pesos) == 0
            pesos = fill(1.0, length(destinos))
        end
    else
        error("modo_prob desconocido: $modo (usa :uniforme_aristas o :ponderado)")
    end
    return (destinos, pesos)
end

"""
    _elegir_companero!(indice_persona, sueltos_por_apellido, destinos, pesos, rng) -> idx_companero | nothing

Dado el conjunto de personas sueltas agrupadas por apellido, elige un apellido
destino segun (destinos, pesos) restringido a apellidos que TODAVIA tengan algun
suelto, y devuelve el indice (en la poblacion) de un companero suelto de ese
apellido, quitandolo de la bolsa de sueltos. Devuelve `nothing` si ningun
destino tiene sueltos disponibles.
"""
function _elegir_companero!(sueltos_por_apellido::Dict{Int,Vector{Int}},
                            destinos::Vector{Int}, pesos::Vector{Float64},
                            rng::AbstractRNG)
    # Filtrar destinos que aun tengan personas sueltas.
    disp_idx = [k for k in eachindex(destinos)
                if haskey(sueltos_por_apellido, destinos[k]) && !isempty(sueltos_por_apellido[destinos[k]])]
    isempty(disp_idx) && return nothing
    dest = destinos[disp_idx]
    w = Weights(pesos[disp_idx])
    apellido_companero = sample(rng, dest, w)
    bolsa = sueltos_por_apellido[apellido_companero]
    # Tomar un suelto cualquiera de ese apellido (el ultimo, barato).
    return pop!(bolsa)
end

"""
    siguiente_generacion_parejas(poblacion, compatibles, rng; modo_prob=:uniforme_aristas)
        -> (nueva_poblacion, parejas_forzadas)

Produce la siguiente generacion segun el modelo por parejas. Devuelve la nueva
poblacion (mismo tamano que la entrada) y un vector `parejas_forzadas` con los
pares de apellidos `(a, b)` que se emparejaron de forma forzada (incompatibles),
es decir, las mutaciones definidas ocurridas en esta generacion.
"""
function siguiente_generacion_parejas(poblacion::Vector{Int},
                                      compatibles::Dict{Int,Vector{Int}},
                                      rng::AbstractRNG = Random.default_rng();
                                      modo_prob::Symbol = :uniforme_aristas)
    N = length(poblacion)
    N == 0 && return (Int[], Tuple{Int,Int}[])
    isodd(N) && error("El modelo por parejas requiere poblacion par (N=$N). Parametriza por numero de parejas.")

    conteos = contar(poblacion)
    vivos = Set(keys(conteos))

    # Pesos de emparejamiento por apellido (calculados una vez por generacion).
    pesos_por_apellido = Dict{Int,Tuple{Vector{Int},Vector{Float64}}}()
    for i in keys(conteos)
        pesos_por_apellido[i] = pesos_emparejamiento(i, compatibles, conteos, vivos, modo_prob)
    end

    # Bolsa de personas sueltas, agrupadas por apellido. Cada persona es un
    # indice 1..N; guardamos su apellido aparte.
    apellido_de = copy(poblacion)
    sueltos_por_apellido = Dict{Int,Vector{Int}}()
    for (idx, a) in enumerate(poblacion)
        push!(get!(sueltos_por_apellido, a, Int[]), idx)
    end

    # Orden aleatorio de procesamiento de las personas.
    orden = shuffle(rng, collect(1:N))

    nueva_poblacion = Int[]
    parejas_forzadas = Tuple{Int,Int}[]

    emparejado = falses(N)

    for persona in orden
        emparejado[persona] && continue
        apellido_i = apellido_de[persona]
        # La persona ya salio de la bolsa en cuanto la procesamos: la quitamos.
        bolsa_i = sueltos_por_apellido[apellido_i]
        pos = findfirst(==(persona), bolsa_i)
        pos !== nothing && deleteat!(bolsa_i, pos)

        destinos, pesos = pesos_por_apellido[apellido_i]
        companero = _elegir_companero!(sueltos_por_apellido, destinos, pesos, rng)

        forzada = false
        if companero === nothing
            # No hay compatible suelto: emparejamiento forzado (mutacion definida).
            companero = _tomar_cualquier_suelto!(sueltos_por_apellido, rng)
            forzada = true
        end

        if companero === nothing
            # No queda nadie mas suelto: esta persona se queda sin pareja.
            # No deberia pasar con N par, pero por robustez la reinsertamos y
            # se emparejara consigo misma (dos hijos de su apellido).
            emparejado[persona] = true
            append!(nueva_poblacion, (apellido_i, apellido_i))
            continue
        end

        emparejado[persona] = true
        emparejado[companero] = true
        apellido_j = apellido_de[companero]

        if forzada
            push!(parejas_forzadas, (apellido_i, apellido_j))
        end

        # Herencia: un miembro al azar pasa su apellido a los DOS hijos.
        apellido_hijos = rand(rng, Bool) ? apellido_i : apellido_j
        append!(nueva_poblacion, (apellido_hijos, apellido_hijos))
    end

    return (nueva_poblacion, parejas_forzadas)
end

"""
    _tomar_cualquier_suelto!(sueltos_por_apellido, rng) -> idx | nothing

Toma una persona suelta cualquiera (sin importar apellido) para un
emparejamiento forzado. Devuelve su indice o `nothing` si no queda nadie.
"""
function _tomar_cualquier_suelto!(sueltos_por_apellido::Dict{Int,Vector{Int}}, rng::AbstractRNG)
    apellidos_con_sueltos = [a for (a, b) in sueltos_por_apellido if !isempty(b)]
    isempty(apellidos_con_sueltos) && return nothing
    a = rand(rng, apellidos_con_sueltos)
    return pop!(sueltos_por_apellido[a])
end

"""
    ResultadoParejas

Resultado de `simular_parejas`:
- `estado_final`: poblacion en la ultima generacion.
- `generaciones`: numero de transiciones simuladas.
- `historial`: Vector con el estado de cada generacion (incluida la 0).
- `forzadas_por_gen`: Vector con la lista de parejas forzadas de cada transicion.
- `absorbido`: `true` si se alcanzo un estado absorbente, ya sea fijacion de un
  unico apellido o un equilibrio permanente de k apellidos (sin parejas mixtas
  posibles entre los vivos).
- `apellidos_final`: apellidos vivos al final.
- `transitorio`: `true` si la corrida se corto por `max_generaciones` y el
  estado todavia admite parejas mixtas (la diversidad aun podria bajar). Es
  excluyente con `absorbido`.
"""
struct ResultadoParejas
    estado_final::Vector{Int}
    generaciones::Int
    historial::Vector{Vector{Int}}
    forzadas_por_gen::Vector{Vector{Tuple{Int,Int}}}
    absorbido::Bool
    apellidos_final::Vector{Int}
    transitorio::Bool
end

"""
    hay_parejas_mixtas_posibles(apellidos_vivos, compatibles) -> Bool

Criterio EXACTO de progreso de la diversidad en el modelo por parejas. La
diversidad solo puede disminuir cuando puede formarse una pareja MIXTA: dos
personas de apellidos DISTINTOS que sean compatibles entre si. Si entre los
apellidos vivos no existe ninguna arista `i -> j` con `i != j` (y reciproca,
porque emparejarse es simetrico), entonces ya no pueden perderse apellidos y el
estado es un equilibrio permanente.

Devuelve `true` si existe al menos un par de apellidos vivos distintos que sean
mutuamente compatibles (hay arista entre ellos en algun sentido).
"""
function hay_parejas_mixtas_posibles(apellidos_vivos, compatibles::Dict{Int,Vector{Int}})
    vivos = Set(apellidos_vivos)
    for i in apellidos_vivos
        vecinos = haskey(compatibles, i) ? compatibles[i] : Int[]
        for j in vecinos
            if j != i && j in vivos
                return true
            end
        end
    end
    return false
end

"""
    simular_parejas(estado_inicial, compatibles; parametros=ParametrosParejas(),
                    rng=Random.default_rng()) -> ResultadoParejas

Itera `siguiente_generacion_parejas` y se detiene cuando ocurre lo primero de:

1. Queda un UNICO apellido (fijacion): `absorbido=true`, `transitorio=false`.
2. Ya no pueden formarse parejas MIXTAS entre los apellidos vivos (criterio
   exacto `hay_parejas_mixtas_posibles`): equilibrio permanente de k apellidos,
   p.ej. Romeo y Julieta. `absorbido=true`, `transitorio=false`.
3. Se alcanza `max_generaciones` sin cumplir 1 ni 2: la corrida se corta. Si a
   esa altura todavia hay parejas mixtas posibles, el estado se reporta como
   `transitorio=true` (la diversidad aun podria bajar con mas generaciones).

Asi, en poblaciones pequenas el proceso corre hasta fijar o hasta un equilibrio
real; en poblaciones grandes, si no se alcanza estabilidad dentro del tope, se
reporta explicitamente como estado transitorio.
"""
function simular_parejas(estado_inicial::Vector{Int},
                         compatibles::Dict{Int,Vector{Int}};
                         parametros::ParametrosParejas = ParametrosParejas(),
                         rng::AbstractRNG = Random.default_rng())
    historial = Vector{Vector{Int}}()
    push!(historial, copy(estado_inicial))
    forzadas_por_gen = Vector{Vector{Tuple{Int,Int}}}()

    poblacion = copy(estado_inicial)
    generaciones = 0
    absorbido = false

    # Estado absorbente inicial (ya no hay parejas mixtas posibles de entrada).
    apellidos_ini = sort(unique(poblacion))
    if length(apellidos_ini) <= 1 || !hay_parejas_mixtas_posibles(apellidos_ini, compatibles)
        return ResultadoParejas(poblacion, 0, historial, forzadas_por_gen,
                                true, apellidos_ini, false)
    end

    for g in 1:parametros.max_generaciones
        poblacion, forzadas = siguiente_generacion_parejas(poblacion, compatibles, rng;
                                                            modo_prob = parametros.modo_prob)
        push!(historial, copy(poblacion))
        push!(forzadas_por_gen, forzadas)
        generaciones += 1

        apellidos_actuales = sort(unique(poblacion))

        # 1. Fijacion: un solo apellido.
        if length(apellidos_actuales) == 1
            absorbido = true
            break
        end

        # 2. Equilibrio permanente: ya no hay parejas mixtas posibles.
        if !hay_parejas_mixtas_posibles(apellidos_actuales, compatibles)
            absorbido = true
            break
        end
    end

    # 3. Si salimos por tope sin absorber, el estado es transitorio si todavia
    #    pueden formarse parejas mixtas (la diversidad aun podria bajar).
    apellidos_final = sort(unique(poblacion))
    transitorio = !absorbido && hay_parejas_mixtas_posibles(apellidos_final, compatibles)

    return ResultadoParejas(poblacion, generaciones, historial, forzadas_por_gen,
                            absorbido, apellidos_final, transitorio)
end
