# scripts/experimentos_parejas.jl
#
# Script REPRODUCIBLE del catalogo de experimentos del modelo por parejas.
# Corre sin Pluto, con semillas fijas, y escribe a disco:
#
#   output/
#     <exp>_grafo.png          grafo inicial de la topologia
#     <exp>_conteos.png        conteo de cada apellido por generacion (una corrida)
#     <exp>_diversidad.png     apellidos vivos por generacion (una corrida)
#     <exp>_supervivientes.png distribucion del nº de supervivientes (muchas corridas)
#     resumen_experimentos.csv resumen agregado de las repeticiones
#
# Ejecutar desde la raiz del proyecto con:
#     julia --project=. scripts/experimentos_parejas.jl
# o:
#     make experimentos
#
# Parametrizable por variables de entorno (igual que run_simulation.jl se corre
# sin argumentos, aqui la configuracion va por el entorno):
#   EXP_N_REPS        nº de repeticiones por experimento        (default 200)
#   EXP_PAREJAS       nº de parejas por apellido en el estado    (default 3)
#   EXP_SEMILLA_BASE  semilla base de las repeticiones           (default 1)
#   EXP_MAX_GEN       tope de generaciones por corrida           (default 1000)
#   EXP_SEMILLA_FIG   semilla de la corrida de ejemplo (figuras) (default 7)
#   EXP_OUTDIR        carpeta de salida                          (default "output")
#
# No agrega dependencias nuevas: usa DelimitedFiles (stdlib) para el CSV.

using Random
using StatsBase
using Printf
using DelimitedFiles

include(joinpath(@__DIR__, "..", "src", "coupleModel.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "graphTopologies.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "experimentMetrics.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "visualization.jl"))

# --------------------------------------------------------------------------
# Lectura de parametros desde el entorno (con defaults reproducibles).
# --------------------------------------------------------------------------

_env_int(clave, def) = haskey(ENV, clave) ? parse(Int, ENV[clave]) : def

struct ConfigExperimentos
    n_reps::Int
    parejas::Int
    semilla_base::Int
    max_gen::Int
    semilla_fig::Int
    outdir::String
end

function leer_config()
    return ConfigExperimentos(
        _env_int("EXP_N_REPS", 200),
        _env_int("EXP_PAREJAS", 3),
        _env_int("EXP_SEMILLA_BASE", 1),
        _env_int("EXP_MAX_GEN", 1000),
        _env_int("EXP_SEMILLA_FIG", 7),
        get(ENV, "EXP_OUTDIR", joinpath(@__DIR__, "..", "output")),
    )
end

# --------------------------------------------------------------------------
# Catalogo fijo de experimentos.
#
# Cada experimento define: un id (prefijo de los archivos), una descripcion, la
# adyacencia `compatibles`, las parejas por apellido de cada apellido y, cuando
# aplica, el apellido `puente` a vigilar. El estado inicial se construye con
# `parejas_por_apellido` x numero de parejas de la config.
# --------------------------------------------------------------------------

struct Experimento
    id::String
    descripcion::String
    compatibles::Dict{Int,Vector{Int}}
    puente::Union{Int,Nothing}
end

"""
    construir_catalogo(parejas) -> Vector{Experimento}

Catalogo fijo de experimentos del modelo por parejas. `parejas` solo se usa para
tamanos que dependan de el; la poblacion concreta se arma luego con
`poblacion_desde_parejas`.
"""
function construir_catalogo(::Int)
    return Experimento[
        Experimento("completo",
                    "Grafo completo de 4 apellidos (baseline: todos compatibles).",
                    topo_completo(4), nothing),
        Experimento("romeo_julieta",
                    "Bipartito aislado {1,2}|{3,4}: dos componentes, equilibrio permanente.",
                    topo_bipartito([1, 2], [3, 4]; cruzado=false), nothing),
        Experimento("cadena_puente",
                    "Cadena 1-2-3-4: los apellidos 2 y 3 son puentes. Vigilamos el puente 2.",
                    topo_cadena(4), 2),
        Experimento("estrella",
                    "Estrella con centro 1 y 4 hojas: el centro es el puente expuesto.",
                    topo_estrella(5; centro=1), 1),
        Experimento("anillo",
                    "Anillo de 4 apellidos: cadena cerrada, sin extremos.",
                    topo_anillo(4), nothing),
        Experimento("grid",
                    "Grid 2x2: vecindad local ortogonal entre 4 apellidos.",
                    topo_grid(2, 2), nothing),
        Experimento("bipartito_cruzado",
                    "Bipartito cruzado {1,2}|{3,4}: cada grupo solo se empareja con el otro.",
                    topo_bipartito([1, 2], [3, 4]; cruzado=true), nothing),
        Experimento("tripartito_ciclico",
                    "Tripartito ciclico {1}-{2}-{3}-{1}: tres apellidos en anillo de grupos.",
                    topo_tripartito([1], [2], [3]; ciclico=true), nothing),
        Experimento("tripartito_cadena",
                    "Tripartito en cadena {1}-{2}-{3}: el apellido 2 es el puente.",
                    topo_tripartito([1], [2], [3]; ciclico=false), 2),
        Experimento("componentes",
                    "Tres componentes aisladas {1,2}|{3,4}|{5,6}: tres linajes permanentes.",
                    topo_componentes([1, 2], [3, 4], [5, 6]), nothing),
    ]
end

estado_inicial_de(exp::Experimento, parejas::Int) =
    poblacion_desde_parejas(Dict(a => parejas for a in sort(collect(keys(exp.compatibles)))))

# --------------------------------------------------------------------------
# Columnas del CSV resumen.
# --------------------------------------------------------------------------

const COLUMNAS_CSV = [
    "experimento",
    "descripcion",
    "n_apellidos",
    "n_reps",
    "parejas_por_apellido",
    "poblacion",
    "modo_prob",
    "max_gen",
    "semilla_base",
    "prob_absorcion",
    "prob_transitorio",
    "gen_media",
    "gen_std",
    "n_vivos_medio",
    "mutaciones_media",
    "prob_puente_primero",
    "dist_n_vivos",
    "prob_ganar",
]

"""
    dict_a_texto(d) -> String

Serializa un Dict{Int,<:Real} a un texto compacto "k:v;k:v" ordenado por clave,
apto para una sola celda de CSV (sin comas, para no romper el delimitador).
"""
function dict_a_texto(d::AbstractDict)
    isempty(d) && return ""
    partes = ["$(k):$(round(float(d[k]); digits=4))" for k in sort(collect(keys(d)))]
    return join(partes, ";")
end

"""
    fila_csv(exp, estado, rep; modo_prob, cfg) -> Vector{Any}

Arma la fila del CSV para un experimento, a partir de la salida agregada de
`correr_repeticiones`.
"""
function fila_csv(exp::Experimento, estado::Vector{Int}, rep; modo_prob::Symbol,
                  cfg::ConfigExperimentos)
    n_vivos_medio = sum(k * v for (k, v) in rep.dist_n_vivos; init = 0.0)
    puente_txt = rep.prob_puente_primero === nothing ? "" :
        string(round(rep.prob_puente_primero; digits = 4))
    return Any[
        exp.id,
        exp.descripcion,
        length(unique(estado)),
        rep.n_reps,
        cfg.parejas,
        length(estado),
        String(modo_prob),
        cfg.max_gen,
        cfg.semilla_base,
        round(rep.prob_absorcion; digits = 4),
        round(rep.prob_transitorio; digits = 4),
        round(rep.gen_media; digits = 4),
        round(rep.gen_std; digits = 4),
        round(n_vivos_medio; digits = 4),
        round(rep.mutaciones_media; digits = 4),
        puente_txt,
        dict_a_texto(rep.dist_n_vivos),
        dict_a_texto(rep.prob_ganar),
    ]
end

# --------------------------------------------------------------------------
# Runner principal.
# --------------------------------------------------------------------------

function main()
    cfg = leer_config()
    modo_prob = :uniforme_aristas

    mkpath(cfg.outdir)

    println("="^70)
    println("Catalogo de experimentos - modelo por parejas")
    println("="^70)
    @printf("  Repeticiones por experimento : %d\n", cfg.n_reps)
    @printf("  Parejas por apellido         : %d\n", cfg.parejas)
    @printf("  Semilla base (repeticiones)  : %d\n", cfg.semilla_base)
    @printf("  Semilla figura (corrida ej.) : %d\n", cfg.semilla_fig)
    @printf("  Max generaciones             : %d\n", cfg.max_gen)
    @printf("  Modo de probabilidad         : %s\n", modo_prob)
    println("  Carpeta de salida            : ", abspath(cfg.outdir))
    println("-"^70)

    catalogo = construir_catalogo(cfg.parejas)
    filas = Vector{Vector{Any}}()

    for exp in catalogo
        estado = estado_inicial_de(exp, cfg.parejas)
        apellidos = sort(collect(keys(exp.compatibles)))
        conteos_ini = contar(estado)

        println("• $(exp.id): $(exp.descripcion)")

        # 1) Corrida de ejemplo (semilla fija) para las figuras de trayectoria.
        resultado = correr_experimento(estado, exp.compatibles;
                                       modo_prob = modo_prob,
                                       max_generaciones = cfg.max_gen,
                                       semilla = cfg.semilla_fig)

        # 2) Repeticiones agregadas (semillas fijas) para el CSV y la distribucion.
        rep = correr_repeticiones(estado, exp.compatibles;
                                  modo_prob = modo_prob,
                                  max_generaciones = cfg.max_gen,
                                  n_reps = cfg.n_reps,
                                  semilla_base = cfg.semilla_base,
                                  puente = exp.puente)

        # 3) Figuras a disco.
        g_grafo = dibujar_grafo(exp.compatibles; conteos = conteos_ini,
                                titulo = "Grafo inicial: $(exp.id)")
        savefig(g_grafo, joinpath(cfg.outdir, "$(exp.id)_grafo.png"))

        g_conteos = graficar_conteos(resultado.historial; apellidos = apellidos,
                                     titulo = "Conteo por apellido: $(exp.id)")
        savefig(g_conteos, joinpath(cfg.outdir, "$(exp.id)_conteos.png"))

        g_div = graficar_diversidad(resultado.historial;
                                    titulo = "Apellidos vivos: $(exp.id)")
        savefig(g_div, joinpath(cfg.outdir, "$(exp.id)_diversidad.png"))

        g_dist = graficar_dist_n_vivos(rep;
                                       titulo = "Supervivientes ($(cfg.n_reps) corridas): $(exp.id)")
        savefig(g_dist, joinpath(cfg.outdir, "$(exp.id)_supervivientes.png"))

        push!(filas, fila_csv(exp, estado, rep; modo_prob = modo_prob, cfg = cfg))

        @printf("    -> absorcion %.0f%%, n_vivos medio %.2f, mutaciones %.2f, %d gen (media)\n",
                100 * rep.prob_absorcion,
                sum(k * v for (k, v) in rep.dist_n_vivos; init = 0.0),
                rep.mutaciones_media, round(Int, rep.gen_media))
    end

    # 4) CSV resumen agregado.
    ruta_csv = joinpath(cfg.outdir, "resumen_experimentos.csv")
    datos = permutedims(hcat(COLUMNAS_CSV))                    # fila de cabeceras
    for f in filas
        datos = vcat(datos, permutedims(hcat(f)))
    end
    writedlm(ruta_csv, datos, ',')

    println("-"^70)
    println("Figuras PNG y CSV escritos en: ", abspath(cfg.outdir))
    println("CSV resumen: ", abspath(ruta_csv))
    println("="^70)
end

main()
