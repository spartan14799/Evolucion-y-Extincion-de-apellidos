# test/test_grafo.jl

using Test
using Printf

include(joinpath(@__DIR__, "..", "src", "graphModel.jl"))
include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))

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

@testset "Validación y Probabilidades Custom (Task 3)" begin
    paso("5/6 Validando grafos correctos e incorrectos ...")
    @testset "Validación de suma y no negatividad" begin
        poblacion = [1, 2]
        grafo = GrafoApellidos(poblacion)
        @test validar_grafo(grafo) == true

        # Inyectar error de suma != 1
        grafo.probabilidades[1][1] = 0.9
        grafo.probabilidades[1][2] = 0.5
        @test_throws ErrorException validar_grafo(grafo)

        # Inyectar error de prob negativa
        grafo.probabilidades[1][1] = -0.2
        grafo.probabilidades[1][2] = 1.2
        @test_throws ErrorException validar_grafo(grafo)
    end

    paso("6/6 Probabilidades asimétricas custom ...")
    @testset "Constructor con probabilidades custom" begin
        poblacion = [1, 2]
        # P(1->1)=0.8, P(1->2)=0.2
        # P(2->1)=0.1, P(2->2)=0.9  (asimétricas)
        probs_custom = Dict(
            1 => Dict(1 => 0.8, 2 => 0.2),
            2 => Dict(1 => 0.1, 2 => 0.9)
        )
        grafo_custom = GrafoApellidos(poblacion; probabilidades=probs_custom)
        
        @test validar_grafo(grafo_custom) == true
        @test probabilidades_salientes(grafo_custom, 1)[2] == 0.2
        @test probabilidades_salientes(grafo_custom, 2)[1] == 0.1
    end
end

# ---------------------------------------------------------------------
# Demo requerida en la Task 3
# ---------------------------------------------------------------------
println("\n--- DEMO TASK 3 ---")
pob_demo3 = [1, 2]
probs_custom = Dict(
    1 => Dict(1 => 0.8, 2 => 0.2),
    2 => Dict(1 => 0.1, 2 => 0.9)
)
println("Construyendo grafo con probabilidades custom asimétricas...")
grafo_valido = GrafoApellidos(pob_demo3; probabilidades=probs_custom)
println("Validando grafo válido: ", validar_grafo(grafo_valido) ? "OK" : "Fallo")

println("\nInyectando un error (suma != 1 en vértice 2)...")
probs_invalidas = Dict(
    1 => Dict(1 => 0.8, 2 => 0.2),
    2 => Dict(1 => 0.1, 2 => 0.5) # Esto suma 0.6
)
grafo_invalido = GrafoApellidos(pob_demo3; probabilidades=probs_invalidas)
try
    validar_grafo(grafo_invalido)
catch e
    println("Error capturado correctamente: ", e.msg)
end
println("-------------------\n")

@testset "Transición de Generación Estocástica (Task 4)" begin
    paso("7/8 Invariantes de la transición ...")
    @testset "Invariantes: tamaño, diversidad, absorbente" begin
        rng = MersenneTwister(42)
        pob_inicial = [1, 1, 2, 2, 3, 3]
        
        pob_nueva = siguiente_generacion_grafo(pob_inicial, rng)
        
        # Tamaño constante (N entra = N sale)
        @test length(pob_nueva) == length(pob_inicial)
        
        # Diversidad no aumenta (los apellidos nuevos son un subconjunto de los antiguos)
        @test issubset(unique(pob_nueva), unique(pob_inicial))
        
        # Estado absorbente (si hay un solo apellido, se mantiene fijo)
        pob_fijada = [2, 2, 2, 2]
        @test siguiente_generacion_grafo(pob_fijada, rng) == [2, 2, 2, 2]
    end
    
    paso("8/8 Reproducibilidad ...")
    @testset "Semilla estocástica (RNG)" begin
        pob = [1, 2, 3, 4, 5]
        res1 = siguiente_generacion_grafo(pob, MersenneTwister(123))
        res2 = siguiente_generacion_grafo(pob, MersenneTwister(123))
        @test res1 == res2
    end
end


# ---------------------------------------------------------------------
# Demo requerida en la Task 4
# ---------------------------------------------------------------------
println("\n--- DEMO TASK 4 ---")
rng_demo = MersenneTwister(100)
pob_t0 = [1, 1, 2, 2, 3, 3]
println("Generación 0: ", pob_t0)
println("Conteos G0:   ", contar_apellidos(pob_t0))

pob_t1 = siguiente_generacion_grafo(pob_t0, rng_demo)
println("Generación 1: ", pob_t1)
println("Conteos G1:   ", contar_apellidos(pob_t1))
println("-------------------\n")

println("\n--- DEMO TASK 5 ---")
rng_demo = MersenneTwister(123)
demo_inicial = [1, 2, 3, 4, 5, 6, 7, 8]
params_demo = ParametrosModelo(max_generaciones=200)

println("Iniciando simulación de fijación para: ", demo_inicial)
resultado = simular_fijacion(demo_inicial; parametros=params_demo, rng=rng_demo)

println("¿Se alcanzó la fijación?: ", resultado.fijacion ? "Sí" : "No")
println("Generaciones transcurridas: ", resultado.generaciones)
println("Apellido ganador: ", resultado.ganador)
println("Historial (length): ", length(resultado.historial))
println("Estado final: ", resultado.historial[end])
println("-------------------\n")
