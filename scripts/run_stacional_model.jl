# scripts/run_stacional_model.jl

using Random
using StatsBase
using Printf

# Incluimos el modelo que ahora usa grafos por debajo
include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))

function main()
    println("="^60)
    println("Simulación de Modelo Estacionario (Grafo Dirigido)")
    println("="^60)

    # 1. Configuración (lee variables de entorno o usa defaults)
    # Ejemplo de uso: POBLACION_INICIAL="1,1,1,2,2,3,3" julia scripts/run_stacional_model.jl
    estado_str = get(ENV, "POBLACION_INICIAL", "1,1,2,2,3,3,4,4")
    estado_inicial = parse.(Int, split(estado_str, ","))
    max_gens = parse(Int, get(ENV, "MAX_GENERACIONES", "100"))

    # Semilla: entero fijo (reproducible) o azar.
    # SEMILLA="nothing" (o variable ausente con default aleatorio) => azar.
    semilla = if haskey(ENV, "SEMILLA")
        ENV["SEMILLA"] == "nothing" ? nothing : parse(Int, ENV["SEMILLA"])
    else
        nothing
    end

    println("Población inicial N = $(length(estado_inicial)): $estado_inicial")
    println("Máximo de generaciones: $max_gens")
    println("Semilla: $(semilla === nothing ? "aleatoria" : semilla)")

    # 2. Imprimir la matriz/tabla de probabilidades de la primera generación
    println("\n--- Grafo de la Generación 0 ---")
    grafo_inicial = GrafoApellidos(estado_inicial)
    println("Vértices y conteos: ", grafo_inicial.conteos)
    println("Tabla de probabilidades salientes P(i -> j):")
    for i in grafo_inicial.apellidos
        print("Desde el apellido $i: ")
        probs = probabilidades_salientes(grafo_inicial, i)
        suma = 0.0
        for (j, p) in sort(collect(probs))
            @printf("%d=%.3f  ", j, p)
            suma += p
        end
        @printf("| Suma: %.1f\n", suma)
    end
    println("--------------------------------\n")

    # 3. Punto de extensión para probabilidades custom (documentado)
    # Para usar probabilidades custom en cada generación, puedes definir una función 
    # y pasarla a ParametrosModelo(fijar_probabilidades=mi_funcion_custom).
    #
    # function mi_estrategia_custom(poblacion_actual::Vector{Int})
    #     # Retornar un Dict{Int, Dict{Int, Float64}} válido
    #     # ...
    # end

    # 4. Configurar y correr la simulación (usando el default proporcional)
    parametros = ParametrosModelo(max_generaciones=max_gens)
    rng = semilla === nothing ? MersenneTwister() : MersenneTwister(semilla)
    
    println("Corriendo simulación estocástica por pareja...")
    resultado = simular_fijacion(estado_inicial; parametros=parametros, rng=rng)

    # 5. Imprimir la evolución usando countmap de StatsBase
    println("\n--- Evolución de los Linajes ---")
    for g in 1:length(resultado.historial)
        # Mostrar las primeras 5 y las últimas 5 generaciones para no saturar la consola
        if g <= 5 || g > length(resultado.historial) - 5
            cmap = countmap(resultado.historial[g])
            str_cmap = join(["$k=>$v" for (k,v) in sort(collect(cmap))], ", ")
            println("Gen $(g-1): $str_cmap")
        elseif g == 6
            println("... (saltando generaciones intermedias) ...")
        end
    end

    # 6. Veredicto final
    println("\n--- Resultado Final ---")
    println("Fijación alcanzada: ", resultado.fijacion ? "Sí" : "No")
    println("Generaciones transcurridas: ", resultado.generaciones)
    if resultado.fijacion
        println("Apellido ganador: ", resultado.ganador)
    end
    println("============================================================")
end

main()
