# test/test_grafo.jl

using Test
using Printf

include(joinpath(@__DIR__, "..", "src", "graphModel.jl"))

const _T0 = Ref(time())
paso(msg) = @printf("[%6.1fs] %s\n", time() - _T0[], msg)

println("="^64)
println("Ejecutando tests de la estructura de Grafo (Task 1)")
println("="^64)
_T0[] = time()

@testset "Estructura del Grafo y Conteo" begin
    
    paso("1/2 Validando conteos y vértices únicos ...")
    @testset "Conteos y vértices" begin
        poblacion = [1, 1, 1, 2, 2, 3, 3]
        grafo = GrafoApellidos(poblacion)
        
        # Validar vértices
        @test grafo.apellidos == [1, 2, 3]
        
        # Validar conteos
        @test grafo.conteos[1] == 3
        @test grafo.conteos[2] == 2
        @test grafo.conteos[3] == 2
    end

    paso("2/2 Validando compatibilidad (grafo completo) ...")
    @testset "Compatibilidad default" begin
        poblacion = [4, 5]
        grafo = GrafoApellidos(poblacion)
        
        @test grafo.compatibles[4] == [4, 5]
        @test grafo.compatibles[5] == [4, 5]
    end
end

# ---------------------------------------------------------------------
# Demo requerida en la Task 1
# ---------------------------------------------------------------------
println("\n--- DEMO TASK 1 ---")
pob_demo = [1, 1, 1, 2, 2, 3, 3]
grafo_demo = GrafoApellidos(pob_demo)

println("Población original: ", pob_demo)
println("Vértices construidos: ", grafo_demo.apellidos)
println("Conteos por vértice:")
letras = Dict(1 => "A", 2 => "B", 3 => "C")
for id in grafo_demo.apellidos
    println("  $(letras[id]) = $(grafo_demo.conteos[id])")
end
println("-------------------\n")

paso("Todos los bloques de Task 1 finalizados.")
