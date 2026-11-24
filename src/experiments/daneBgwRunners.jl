# src/experiments/daneBgwRunners.jl
#
# Orquestacion Monte Carlo del MODELO BGW aplicada a los datos del Censo DANE
# 2018. Alimenta el notebook notebooks/results/bgw_colombia_dane.jl: separa la
# logica de simulacion de la presentacion, de modo que el notebook solo dibuja.
# Llama a `simular_generaciones` de src/bgwModel.jl; no reimplementa el modelo.
#
# Documentacion completa (funciones, uso): docs/code/daneBgwRunners.md
#
# Requiere cargados antes: src/core/numericalMethods.jl, src/bgwModel.jl y
# src/experiments/daneData.jl.

using Random
using StatsBase

"""
    correr_montecarlo_colombia(modelo, z0, n_sims, gen_eval, semilla) -> NamedTuple

Monte Carlo del tiempo de extincion de un apellido bajo el modelo BGW, partiendo
de `z0` portadores (p. ej. el conteo DANE de un apellido). Devuelve un NamedTuple:
- `tiempos_ext`          : generacion de extincion de cada corrida que se extinguio.
- `trayectorias_muestra` : hasta 6 trayectorias completas para graficar.
- `extintos_en_X`        : cuantas corridas se extinguieron en <= `gen_eval` generaciones.
- `prob_ext_X`           : fraccion de corridas extintas en <= `gen_eval`.
- `mediana_ext`, `media_ext` : estadisticos del tiempo de extincion.
- `max_gen_sim`          : horizonte de generaciones efectivamente simulado.
"""
function correr_montecarlo_colombia(modelo::ModeloPBGW, z0::Int, n_sims::Int,
                                    gen_eval::Int, semilla::Int)
    rng = MersenneTwister(semilla)
    max_gen_sim = max(gen_eval + 15, 35)
    tiempos_ext = Int[]
    trayectorias_muestra = Vector{Vector{Int}}()
    extintos_en_X = 0

    for s in 1:n_sims
        traj = simular_generaciones(modelo, max_gen_sim; poblacion_inicial = z0, rng = rng)
        if s <= 6
            push!(trayectorias_muestra, copy(traj))
        end
        idx_cero = findfirst(==(0), traj)
        if idx_cero !== nothing
            gen_ext = idx_cero - 1
            push!(tiempos_ext, gen_ext)
            if gen_ext <= gen_eval
                extintos_en_X += 1
            end
        end
    end

    prob_ext_X = extintos_en_X / n_sims
    mediana_ext = isempty(tiempos_ext) ? 0 : round(Int, median(tiempos_ext))
    media_ext = isempty(tiempos_ext) ? 0.0 : mean(tiempos_ext)

    return (tiempos_ext = tiempos_ext,
            trayectorias_muestra = trayectorias_muestra,
            extintos_en_X = extintos_en_X,
            prob_ext_X = prob_ext_X,
            mediana_ext = mediana_ext,
            media_ext = media_ext,
            max_gen_sim = max_gen_sim)
end

"""
    calcular_medianas_10(modelo; n_reps=30, semilla=123) -> Vector{Float64}

Para cada uno de los 10 apellidos del Censo DANE (poblacion inicial =
`CONTEOS_DANE[i]`), estima por Monte Carlo la generacion MEDIANA de extincion
bajo `modelo`, con `n_reps` corridas por apellido. Devuelve el vector de medianas
en el mismo orden que `APELLIDOS_DANE`.
"""
function calcular_medianas_10(modelo::ModeloPBGW; n_reps::Int = 30, semilla::Int = 123)
    tiempos_medios = Float64[]
    rng_comp = MersenneTwister(semilla)
    for cnt in CONTEOS_DANE
        tiempos_i = Int[]
        for _ in 1:n_reps
            traj_i = simular_generaciones(modelo, 35; poblacion_inicial = cnt, rng = rng_comp)
            idx_0 = findfirst(==(0), traj_i)
            push!(tiempos_i, idx_0 !== nothing ? idx_0 - 1 : 36)
        end
        push!(tiempos_medios, median(tiempos_i))
    end
    return tiempos_medios
end
