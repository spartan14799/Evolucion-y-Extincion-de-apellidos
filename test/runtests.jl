# test/runtests.jl
#
# Suite de pruebas de correctitud para los metodos numericos.
#
# Ejecutar desde la raiz del proyecto con:
#     julia test/runtests.jl
#
# Usa unicamente la libreria estandar `Test`, por lo que no requiere
# instalar dependencias externas.

using Test

# Carga las funciones a probar (biseccion, newton, secante).
include(joinpath(@__DIR__, "..", "src", "numericalMethods.jl"))

# Tolerancia con la que comparamos la raiz encontrada contra la raiz exacta.
# Es un orden de magnitud mas holgada que el `epsilon`/`delta` por defecto
# de los metodos (1e-5) para absorber el error de aproximacion.
const TOL = 1e-4

@testset "Metodos numericos" begin

    # ---------------------------------------------------------------
    # Funciones de prueba con raices conocidas analiticamente.
    # ---------------------------------------------------------------
    # f1(x) = x^2 - 2      -> raiz en sqrt(2) ≈ 1.4142135
    f1(x)  = x^2 - 2
    df1(x) = 2x

    # f2(x) = cos(x) - x    -> raiz en ≈ 0.7390851
    f2(x)  = cos(x) - x
    df2(x) = -sin(x) - 1

    # f3(x) = x^3 - x - 2   -> raiz en ≈ 1.5213797
    f3(x)  = x^3 - x - 2
    df3(x) = 3x^2 - 1

    raiz_f1 = sqrt(2)
    raiz_f2 = 0.7390851332151607
    raiz_f3 = 1.5213797068045676

    @testset "Biseccion" begin
        r1 = biseccion(f1, 0.0, 2.0)
        @test r1.convergio
        @test isapprox(r1.raiz, raiz_f1; atol = TOL)
        @test f1(r1.raiz) ≈ 0.0 atol = 1e-3

        r2 = biseccion(f2, 0.0, 1.0)
        @test r2.convergio
        @test isapprox(r2.raiz, raiz_f2; atol = TOL)

        r3 = biseccion(f3, 1.0, 2.0)
        @test r3.convergio
        @test isapprox(r3.raiz, raiz_f3; atol = TOL)

        # Debe rechazar intervalos sin cambio de signo.
        @test_throws ErrorException biseccion(f1, 2.0, 3.0)
    end

    @testset "Newton" begin
        r1 = newton(f1, df1, 1.0)
        @test r1.convergio
        @test isapprox(r1.raiz, raiz_f1; atol = TOL)

        r2 = newton(f2, df2, 0.5)
        @test r2.convergio
        @test isapprox(r2.raiz, raiz_f2; atol = TOL)

        r3 = newton(f3, df3, 1.5)
        @test r3.convergio
        @test isapprox(r3.raiz, raiz_f3; atol = TOL)

        # Si x0 ya es (casi) la raiz, debe converger en 0 iteraciones.
        r0 = newton(f1, df1, raiz_f1)
        @test r0.convergio
        @test r0.iteraciones == 0

        # Derivada nula debe lanzar error (g'(0) = 0 para g(x)=x^2-2 ... 2*0=0).
        @test_throws ErrorException newton(f1, df1, 0.0)
    end

    @testset "Secante" begin
        r1 = secante(f1, 0.0, 2.0)
        @test r1.convergio
        @test isapprox(r1.raiz, raiz_f1; atol = TOL)

        r2 = secante(f2, 0.0, 1.0)
        @test r2.convergio
        @test isapprox(r2.raiz, raiz_f2; atol = TOL)

        r3 = secante(f3, 1.0, 2.0)
        @test r3.convergio
        @test isapprox(r3.raiz, raiz_f3; atol = TOL)
    end

    @testset "Consistencia entre metodos" begin
        # Los tres metodos deben converger a la misma raiz para f1.
        rb = biseccion(f1, 0.0, 2.0)
        rn = newton(f1, df1, 1.0)
        rs = secante(f1, 0.0, 2.0)
        @test isapprox(rb.raiz, rn.raiz; atol = TOL)
        @test isapprox(rn.raiz, rs.raiz; atol = TOL)
    end
end
