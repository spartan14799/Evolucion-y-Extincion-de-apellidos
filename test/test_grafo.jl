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

@testset "Cálculo de probabilidades (Task 2)" begin
    paso("3/4 Calculando probabilidades (grafo completo) ...")
    @testset "Probabilidades proporcionales default" begin
        # A=3, B=2, C=2 => Total 7
        poblacion = [1, 1, 1, 2, 2, 3, 3]
        grafo = GrafoApellidos(poblacion)
        
        # Inyectamos el resultado en el grafo
        grafo.probabilidades = probabilidades_proporcionales(grafo)
        
        probs_A = probabilidades_salientes(grafo, 1)
        @test isapprox(probs_A[1], 3/7)
        @test isapprox(probs_A[2], 2/7)
        @test isapprox(probs_A[3], 2/7)
        
        # Verificar suma = 1 por nodo
        for i in grafo.apellidos
            suma_salientes = sum(values(probabilidades_salientes(grafo, i)))
            @test isapprox(suma_salientes, 1.0)
        end
    end

    paso("4/4 Probabilidades con compatibles restringidos ...")
    @testset "Probabilidades con restricciones" begin
        poblacion = [1, 1, 1, 2, 2, 3, 3]
        
        # A solo compatible con A y B (3 + 2 = 5 personas)
        compatibles_restringidos = Dict(
            1 => [1, 2],       
            2 => [2, 3],       
            3 => [1, 2, 3]     
        )
        grafo = GrafoApellidos(poblacion; compatibles=compatibles_restringidos)
        grafo.probabilidades = probabilidades_proporcionales(grafo)
        
        probs_A = probabilidades_salientes(grafo, 1)
        @test isapprox(probs_A[1], 3/5)
        @test isapprox(probs_A[2], 2/5)
        @test !haskey(probs_A, 3) # No debe existir arista A -> C
        @test isapprox(sum(values(probs_A)), 1.0)
    end
end

# ---------------------------------------------------------------------
# Demo requerida en la Task 2
# ---------------------------------------------------------------------
println("\n--- DEMO TASK 2 ---")
grafo_demo.probabilidades = probabilidades_proporcionales(grafo_demo)
println("Tabla de probabilidades (P i -> j):")
for i in grafo_demo.apellidos
    print("Desde $(letras[i]): ")
    probs = probabilidades_salientes(grafo_demo, i)
    suma = 0.0
    for (j, p) in sort(collect(probs)) # sort para mantener el orden de impresión
        @printf("%s=%.3f  ", letras[j], p)
        suma += p
    end
    @printf("| Suma salientes: %.1f\n", suma)
end
println("-------------------\n")
