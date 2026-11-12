# scripts/run_stacional_model.jl
#
# Ejecuta una simulacion del modelo estacionario (poblacion constante) y muestra
# la evolucion de los apellidos generacion por generacion hasta la fijacion.
#
# Ejecutar desde la raiz del proyecto con:
#     julia --project=. scripts/run_stacional_model.jl
#
# -----------------------------------------------------------------------------
# COMO CAMBIAR LOS PARAMETROS
# -----------------------------------------------------------------------------
# Edita el bloque "PARAMETROS DE LA SIMULACION" dentro de `main()` mas abajo, o
# exporta variables de entorno antes de ejecutar (tienen prioridad sobre los
# valores por defecto del codigo):
#
#     POBLACION_INICIAL=20 APELLIDOS_INICIALES=6 SEMILLA=7 \
#         julia --project=. scripts/run_stacional_model.jl
#
# Variables reconocidas:
#   POBLACION_INICIAL    Numero de individuos (debe ser par). Def: 10
#   APELLIDOS_INICIALES  Apellidos distintos iniciales (<= poblacion). Def: 4
#   SEMILLA              Entero para reproducibilidad, o "nothing" para azar. Def: 42
#   TAMANO_PAREJA        Individuos por grupo que eligen un apellido. Def: 2
#   MAX_GENERACIONES     Tope de seguridad de generaciones. Def: 10000

using Random
using StatsBase

include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))

"""
    configurar_poblacion(tamano, num_apellidos; tamano_pareja = 2, rng)

Construye una poblacion inicial de `tamano` individuos con `num_apellidos`
apellidos distintos (ids `1:num_apellidos`), rellenando el resto al azar.

Valida que `tamano` sea divisible por `tamano_pareja` y que no haya mas
apellidos que individuos.
"""
function configurar_poblacion(tamano::Int, num_apellidos::Int;
                              tamano_pareja::Int = 2,
                              rng = Random.default_rng())
    if tamano % tamano_pareja != 0
        error("La poblacion inicial ($tamano) debe ser divisible por el tamano de pareja ($tamano_pareja).")
    end
    if num_apellidos > tamano
        error("No puede haber mas apellidos ($num_apellidos) que personas ($tamano).")
    end
    if num_apellidos < 1
        error("Debe haber al menos un apellido.")
    end
    poblacion = collect(1:num_apellidos)
    faltantes = tamano - num_apellidos
    if faltantes > 0
        append!(poblacion, rand(rng, 1:num_apellidos, faltantes))
    end

    return shuffle(rng, poblacion)
end

# Lee una variable de entorno como Int, con valor por defecto.
_env_int(nombre, def) = haskey(ENV, nombre) ? parse(Int, ENV[nombre]) : def

function main()
    # ==========================================
    # PARAMETROS DE LA SIMULACION
    # (editables aqui, o via variables de entorno)
    # ==========================================
    POBLACION_INICIAL   = _env_int("POBLACION_INICIAL", 10)    # Debe ser par
    APELLIDOS_INICIALES = _env_int("APELLIDOS_INICIALES", 4)   # <= poblacion
    TAMANO_PAREJA       = _env_int("TAMANO_PAREJA", 2)
    MAX_GENERACIONES    = _env_int("MAX_GENERACIONES", 10_000)

    # Semilla: entero fijo (reproducible) o `nothing` (azar).
    SEMILLA = if haskey(ENV, "SEMILLA")
        ENV["SEMILLA"] == "nothing" ? nothing : parse(Int, ENV["SEMILLA"])
    else
        42
    end

    rng = SEMILLA !== nothing ? MersenneTwister(SEMILLA) : Random.default_rng()
    parametros = ParametrosModelo(tamano_pareja = TAMANO_PAREJA,
                                  max_generaciones = MAX_GENERACIONES)

    println("="^64)
    println("SIMULACION DE MARKOV: POBLACION CONSTANTE")
    println("Individuos: $POBLACION_INICIAL | Apellidos iniciales distintos: $APELLIDOS_INICIALES")
    println("Tamano de pareja: $TAMANO_PAREJA | Semilla: $(SEMILLA === nothing ? "aleatoria" : SEMILLA)")
    println("="^64)

    # Configurar estado inicial.
    estado_inicial = configurar_poblacion(POBLACION_INICIAL, APELLIDOS_INICIALES;
                                           tamano_pareja = TAMANO_PAREJA, rng = rng)

    # Ejecutar la simulacion.
    resultado = simular_fijacion(estado_inicial; parametros = parametros, rng = rng)

    # Mostrar la evolucion estado por estado.
    println("\n--- EVOLUCION POR GENERACION ---")
    for (gen, estado) in enumerate(resultado.historial)
        n_gen = gen - 1
        distribucion = countmap(estado)

        resumen = join(["ID $k:$v" for (k, v) in sort(collect(distribucion))], ", ")
        println("Gen $n_gen | Estado: $estado")
        println("        ↳ Frecuencias: [$resumen]")
    end

    println("-"^64)
    if resultado.fijacion
        println("RESULTADO FINAL: FIJACION ALCANZADA")
        println("El apellido dominante fue el ID: $(resultado.ganador)")
        println("Generaciones hasta la extincion del resto: $(resultado.generaciones)")
    else
        println("RESULTADO FINAL: LIMITE ALCANZADO")
        println("La simulacion se detuvo por seguridad tras $(resultado.generaciones) generaciones.")
        println("Aun sobrevivian varios apellidos.")
    end
    println("="^64)
end

main()
