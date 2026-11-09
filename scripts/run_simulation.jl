# scripts/run_simulation.jl
#
# Simulacion inicial de la evolucion y extincion de un apellido mediante
# un proceso de Galton-Watson (ModeloPBGW).
#
# Ejecutar desde la raiz del proyecto con:
#     julia scripts/run_simulation.jl

using Random
using Printf

include(joinpath(@__DIR__, "..", "src", "numericalMethods.jl"))
include(joinpath(@__DIR__, "..", "src", "markovModels.jl"))

function resumen_modelo(nombre, modelo)
    m = esperanza(modelo)
    q = probabilidad_extincion(modelo)
    println("-"^60)
    println("Modelo: $nombre")
    @printf("  Esperanza (hijos por individuo) m = %.4f\n", m)
    println("  Regimen: ", m < 1 ? "subcritico" : (m == 1 ? "critico" : "supercritico"))
    @printf("  Probabilidad de extincion teorica q = %.4f\n", q)
    return q
end

function main()
    rng = MersenneTwister()  # semilla fija para reproducibilidad

    println("="^60)
    println("Evolucion y extincion de apellidos - Simulacion inicial")
    println("="^60)

    # Modelo supercritico: natalidad alta (m > 1), extincion posible pero < 1.
    modelo_alto = ModeloPBGW(1.5)
    q_alto = resumen_modelo("Natalidad alta (Poisson lambda=1.5)", modelo_alto)

    # Modelo critico/subcritico: natalidad baja (m <= 1), extincion segura.
    modelo_bajo = ModeloPBGW(0.8)
    resumen_modelo("Natalidad baja (Poisson lambda=0.8)", modelo_bajo)

    # Trayectorias de ejemplo para el modelo supercritico.
    println("-"^60)
    println("Trayectorias de ejemplo (modelo natalidad alta, 15 generaciones):")
    for i in 1:5
        traj = simular_generaciones(modelo_alto, 15; poblacion_inicial=1, rng=rng)
        estado = traj[end] == 0 ? "EXTINTO" : "vivo ($(traj[end]))"
        println("  Sim $i: ", traj, "  -> ", estado)
    end

    # Validacion Monte Carlo contra la probabilidad teorica.
    # Usamos un tope de poblacion para que los linajes supercriticos que
    # explotan (crecimiento exponencial) no vuelvan la simulacion inviable.
    println("-"^60)
    println("Validacion Monte Carlo (modelo natalidad alta):")
    gens  = 40
    n_sim = 10_000
    freq = frecuencia_extincion(modelo_alto, gens;
                                n_sim=n_sim, tope_poblacion=10_000, rng=rng)
    @printf("  Frecuencia de extincion simulada (%d gen, n=%d) = %.4f\n", gens, n_sim, freq)
    @printf("  Probabilidad de extincion teorica             = %.4f\n", q_alto)
    @printf("  Diferencia                                    = %.4f\n", abs(freq - q_alto))
    println("="^60)
end

main()
