# test/test_stacional.jl
#
# Pruebas de correctitud para el modelo estacionario (poblacion constante):
# `siguiente_generacion` y `simular_fijacion` de src/stacional-model.jl.
#
# Ejecutar desde la raiz del proyecto con:
#     julia --project=. test/test_stacional.jl
#
# Usa unicamente la libreria estandar `Test` y `Random`.

using Test
using Random
using Printf

include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))

# Utilidad para mostrar el progreso en vivo (mismo estilo que test_markov.jl).
const _T0 = Ref(time())
paso(msg) = @printf("[%6.1fs] %s\n", time() - _T0[], msg)

println("="^64)
println("Ejecutando tests del modelo estacionario (poblacion constante)")
println("="^64)
_T0[] = time()

@testset "Modelo estacionario (poblacion constante)" begin

    paso("1/6  Parametros del modelo ...")
    @testset "ParametrosModelo" begin
        # Valores por defecto.
        p = ParametrosModelo()
        @test p.tamano_pareja == 2
        @test p.max_generaciones == 10_000

        # Constructor por palabras clave.
        p2 = ParametrosModelo(tamano_pareja = 3, max_generaciones = 500)
        @test p2.tamano_pareja == 3
        @test p2.max_generaciones == 500
    end

    paso("2/6  siguiente_generacion: conserva el tamano ...")
    @testset "siguiente_generacion conserva el tamano" begin
        rng = MersenneTwister(1)
        pobl = collect(1:10)                      # tamano par, divisible por 2
        nueva = siguiente_generacion(pobl, rng)
        @test length(nueva) == length(pobl)       # INVARIANTE: tamano constante

        # La diversidad nunca aumenta: todo apellido nuevo ya existia.
        @test issubset(Set(nueva), Set(pobl))

        # Cada pareja (bloque de 2) comparte el mismo apellido.
        for i in 1:2:length(nueva)
            @test nueva[i] == nueva[i + 1]
        end
    end

    paso("3/6  siguiente_generacion: tamano de pareja configurable ...")
    @testset "siguiente_generacion con tamano_pareja = 4" begin
        rng = MersenneTwister(2)
        params = ParametrosModelo(tamano_pareja = 4)
        pobl = collect(1:12)                      # divisible por 4
        nueva = siguiente_generacion(pobl, rng; parametros = params)
        @test length(nueva) == 12

        # Cada bloque de 4 comparte el mismo apellido.
        for i in 1:4:length(nueva)
            @test length(unique(nueva[i:i+3])) == 1
        end

        # Poblacion no divisible por el tamano de pareja => error.
        @test_throws ErrorException siguiente_generacion(collect(1:10), rng; parametros = params)

        # Tamano de pareja invalido => error.
        @test_throws ErrorException siguiente_generacion(collect(1:4), rng;
                                        parametros = ParametrosModelo(tamano_pareja = 0))
    end

    paso("4/6  siguiente_generacion: reproducibilidad ...")
    @testset "Reproducibilidad con misma semilla" begin
        pobl = collect(1:20)
        a = siguiente_generacion(pobl, MersenneTwister(99))
        b = siguiente_generacion(pobl, MersenneTwister(99))
        @test a == b
    end

    paso("5/6  simular_fijacion: estado absorbente y limites ...")
    @testset "simular_fijacion" begin
        # Poblacion ya fijada: 0 generaciones, ganador inmediato.
        res0 = simular_fijacion([7, 7, 7, 7]; rng = MersenneTwister(3))
        @test res0.fijacion
        @test res0.ganador == 7
        @test res0.generaciones == 0
        @test length(res0.historial) == 1

        # El historial siempre empieza con el estado inicial.
        inicial = [1, 2, 3, 4]
        res = simular_fijacion(inicial; rng = MersenneTwister(5))
        @test res.historial[1] == inicial

        # Todos los estados del historial conservan el tamano de la poblacion.
        @test all(e -> length(e) == length(inicial), res.historial)

        # El historial tiene generaciones + 1 estados.
        @test length(res.historial) == res.generaciones + 1

        # Si hay fijacion, el ganador aparece en todo el estado final.
        if res.fijacion
            @test all(==(res.ganador), res.historial[end])
            @test length(unique(res.historial[end])) == 1
        end

        # Limite de generaciones: con max_generaciones = 0 no avanza.
        res_lim = simular_fijacion([1, 2, 3, 4];
                                   max_generaciones = 0, rng = MersenneTwister(1))
        @test res_lim.generaciones == 0
        @test !res_lim.fijacion
        @test res_lim.ganador === nothing
    end

    paso("6/6  simular_fijacion: siempre converge a fijacion ...")
    @testset "Convergencia a fijacion (poblacion pequena)" begin
        # Para poblaciones pequenas la fijacion es casi segura; con un tope
        # amplio debe alcanzarse siempre. Probamos varias semillas.
        for s in 1:20
            res = simular_fijacion([1, 2, 3, 4, 5, 6, 7, 8];
                                   max_generaciones = 10_000,
                                   rng = MersenneTwister(s))
            @test res.fijacion
            @test res.ganador in 1:8
        end
    end
end

paso("Todos los bloques finalizados.")
