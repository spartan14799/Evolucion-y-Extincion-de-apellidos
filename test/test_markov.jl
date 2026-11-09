# test/test_markov.jl
#
# Pruebas de correctitud para el modelo de Galton-Watson (ModeloPBGW) y la
# simulacion de evolucion/extincion de apellidos.
#
# Ejecutar desde la raiz del proyecto con:
#     julia --project=. test/test_markov.jl
#
# Requiere el paquete `Distributions` (ya usado por src/markovModels.jl).

using Test
using Random
using Printf

include(joinpath(@__DIR__, "..", "src", "numericalMethods.jl"))
include(joinpath(@__DIR__, "..", "src", "markovModels.jl"))

# ---------------------------------------------------------------------
# Utilidad para mostrar el progreso en vivo. Imprime el nombre del bloque
# y el tiempo que tardo, para que se vea que el test avanza y no esta
# colgado (util cuando el Monte Carlo tarda unos segundos).
# ---------------------------------------------------------------------
const _T0 = Ref(time())
paso(msg) = @printf("[%6.1fs] %s\n", time() - _T0[], msg)

println("="^64)
println("Ejecutando tests del modelo de Markov / Galton-Watson")
println("="^64)
_T0[] = time()

@testset "Modelo de Markov / Galton-Watson" begin

    paso("1/7  Construccion de ModeloPBGW ...")
    @testset "Construccion de ModeloPBGW" begin
        # Constructor por vector de probabilidades.
        m = ModeloPBGW([0.25, 0.25, 0.5])
        @test length(m.pk) == 3
        @test isapprox(sum(m.pk), 1.0; atol = 1e-9)

        # Rechaza probabilidades negativas.
        @test_throws ErrorException ModeloPBGW([-0.1, 0.6, 0.5])

        # Rechaza vectores que no suman 1.
        @test_throws ErrorException ModeloPBGW([0.2, 0.2, 0.2])

        # Constructor por tasa de natalidad (Poisson) debe normalizar a 1.
        mp = ModeloPBGW(1.5)
        @test isapprox(sum(mp.pk), 1.0; atol = 1e-9)
    end

    paso("2/7  Esperanza (numero medio de hijos) ...")
    @testset "Esperanza (numero medio de hijos)" begin
        # pk = [P(0), P(1), P(2)] = [0.25, 0.25, 0.5]
        # m = 0*0.25 + 1*0.25 + 2*0.5 = 1.25
        m = ModeloPBGW([0.25, 0.25, 0.5])
        @test isapprox(esperanza(m), 1.25; atol = 1e-9)

        # Para Poisson(lambda), m ~ lambda (tras truncar y normalizar).
        mp = ModeloPBGW(2.0; max_hijos = 30)
        @test isapprox(esperanza(mp), 2.0; atol = 1e-3)
    end

    paso("3/7  Criterio de extincion segura (m <= 1) ...")
    @testset "Criterio de extincion segura (m <= 1)" begin
        subcritico = ModeloPBGW([0.6, 0.3, 0.1])   # m = 0.5
        critico    = ModeloPBGW([0.5, 0.0, 0.5])   # m = 1.0
        supercrit  = ModeloPBGW([0.2, 0.3, 0.5])   # m = 1.3

        @test es_extincion_segura(subcritico)
        @test es_extincion_segura(critico)
        @test !es_extincion_segura(supercrit)
    end

    paso("4/7  Funcion generadora de probabilidad (fgp) ...")
    @testset "Funcion generadora de probabilidad (fgp)" begin
        m = ModeloPBGW([0.25, 0.25, 0.5])
        # g(1) = suma de probabilidades = 1.
        @test isapprox(fgp(m, 1.0), 1.0; atol = 1e-9)
        # g(0) = P(0 hijos) = pk[1].
        @test isapprox(fgp(m, 0.0), m.pk[1]; atol = 1e-9)
        # g(t) = 0.25 + 0.25 t + 0.5 t^2  -> en t=0.5: 0.25+0.125+0.125 = 0.5
        @test isapprox(fgp(m, 0.5), 0.5; atol = 1e-9)
    end

    paso("5/7  Probabilidad de extincion q (biseccion y secante) ...")
    @testset "Probabilidad de extincion q" begin
        # Caso subcritico/critico: q = 1.
        @test probabilidad_extincion(ModeloPBGW([0.6, 0.3, 0.1])) == 1.0
        @test probabilidad_extincion(ModeloPBGW([0.5, 0.0, 0.5])) == 1.0

        # Caso supercritico con solucion analitica.
        # pk = [p0, 0, p2] con p2 = 0.6, p0 = 0.4, m = 1.2 > 1.
        # g(t) = 0.4 + 0.6 t^2 = t  =>  0.6 t^2 - t + 0.4 = 0
        #   t = (1 ± sqrt(1 - 0.96)) / 1.2 = (1 ± 0.2)/1.2
        #   menor raiz = 0.8/1.2 = 0.6666...
        m = ModeloPBGW([0.4, 0.0, 0.6])
        q = probabilidad_extincion(m)
        @test isapprox(q, 2/3; atol = 1e-4)
        @test 0.0 <= q <= 1.0

        # El metodo de la secante debe dar la misma q.
        q_sec = probabilidad_extincion(m; buscador = :secante)
        @test isapprox(q_sec, 2/3; atol = 1e-3)

        # Buscador invalido lanza error.
        @test_throws ErrorException probabilidad_extincion(m; buscador = :foo)
    end

    paso("6/7  Simulacion de generaciones ...")
    @testset "Simulacion de generaciones" begin
        rng = MersenneTwister(123)
        m = ModeloPBGW(1.5)

        traj = simular_generaciones(m, 10; poblacion_inicial = 1, rng = rng)
        @test length(traj) == 11               # generaciones + 1
        @test traj[1] == 1                      # poblacion inicial
        @test all(x -> x >= 0, traj)            # poblacion nunca negativa

        # Una vez en 0, permanece en 0 (absorbente).
        idx = findfirst(==(0), traj)
        if idx !== nothing
            @test all(traj[idx:end] .== 0)
        end

        # Poblacion inicial 0 => extinto en toda la trayectoria.
        traj0 = simular_generaciones(m, 5; poblacion_inicial = 0, rng = rng)
        @test all(traj0 .== 0)

        # 0 generaciones => solo la poblacion inicial.
        traj_cero = simular_generaciones(m, 0; poblacion_inicial = 3, rng = rng)
        @test traj_cero == [3]

        # Argumentos invalidos.
        @test_throws ErrorException simular_generaciones(m, -1)
        @test_throws ErrorException simular_generaciones(m, 5; poblacion_inicial = -2)

        # Reproducibilidad con misma semilla.
        a = simular_generaciones(m, 20; rng = MersenneTwister(7))
        b = simular_generaciones(m, 20; rng = MersenneTwister(7))
        @test a == b
    end

    paso("7/7  Monte Carlo vs probabilidad teorica (puede tardar unos seg) ...")
    @testset "Monte Carlo vs probabilidad teorica" begin
        rng = MersenneTwister(2024)
        m = ModeloPBGW([0.4, 0.0, 0.6])        # q teorica = 2/3
        q_teo = probabilidad_extincion(m)

        # Suficientes simulaciones/generaciones para acercarse a la extincion
        # asintotica, manteniendo el test rapido (~1-2 s).
        n_sim = 8_000
        gens  = 60
        freq = frecuencia_extincion(m, gens; n_sim = n_sim, rng = rng)
        paso(@sprintf("      q_teorica = %.4f  |  freq_simulada = %.4f  (n=%d, gen=%d)",
                      q_teo, freq, n_sim, gens))
        @test isapprox(freq, q_teo; atol = 0.03)  # margen estadistico
    end
end

paso("Todos los bloques finalizados.")
