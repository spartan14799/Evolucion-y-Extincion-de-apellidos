# scripts/run_stacional_model.jl

using Random
using StatsBase

include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))

function configurar_poblacion(tamano::Int, num_apellidos::Int; rng=Random.default_rng())
    if tamano % 2 != 0
        error("La población inicial debe ser un número par para formar parejas exactas.")
    end
    if num_apellidos > tamano
        error("No puede haber más apellidos que personas.")
    end
    poblacion = collect(1:num_apellidos)
    faltantes = tamano - num_apellidos
    if faltantes > 0
        append!(poblacion, rand(rng, 1:num_apellidos, faltantes))
    end
    
    return shuffle(rng, poblacion)
end

function main()
    # ==========================================
    # PARÁMETROS DE LA SIMULACIÓN
    # ==========================================
    POBLACION_INICIAL = 10       # Debe ser par
    APELLIDOS_INICIALES = 4      # Máximo igual a la población
    SEMILLA = 42                 # Cambia el número para resultados distintos, o usa 'nothing'
    
    rng = SEMILLA !== nothing ? MersenneTwister(SEMILLA) : Random.default_rng()

    println("="^64)
    println("SIMULACIÓN DE MARKOV: POBLACIÓN CONSTANTE")
    println("Individuos: $POBLACION_INICIAL | Apellidos iniciales distintos: $APELLIDOS_INICIALES")
    println("="^64)

    # Configurar estado inicial
    estado_inicial = configurar_poblacion(POBLACION_INICIAL, APELLIDOS_INICIALES; rng=rng)
    
    # Ejecutar la simulación
    resultado = simular_fijacion(estado_inicial; rng=rng)

    # Mostrar la evolución estado por estado
    println("\n--- EVOLUCIÓN POR GENERACIÓN ---")
    for (gen, estado) in enumerate(resultado.historial)
        n_gen = gen - 1
        distribucion = countmap(estado)
        
        resumen = join(["ID $k:$v" for (k, v) in sort(collect(distribucion))], ", ")
        println("Gen $n_gen | Estado: $estado")
        println("        ↳ Frecuencias: [$resumen]")
    end

    # ... (código anterior) ...

    println("-"^64)
    if resultado.fijacion
        println("RESULTADO FINAL: FIJACIÓN ALCANZADA")
        println("El apellido dominante fue el ID: $(resultado.ganador)")
        println("Generaciones hasta la extinción del resto: $(resultado.generaciones)")
    else
        println("RESULTADO FINAL: LÍMITE ALCANZADO")
        println("La simulación se detuvo por seguridad tras $(resultado.generaciones) generaciones.")
        println("Aún sobrevivían varios apellidos.")
    end
    println("="^64)

end

main()
