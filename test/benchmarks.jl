# test/benchmarks.jl
#
# Benchmark de los metodos numericos (biseccion, newton, secante).
#
# Ejecutar desde la raiz del proyecto con:
#     julia test/benchmarks.jl
#
# Si el paquete `BenchmarkTools` esta instalado, se usa para obtener
# estadisticas robustas (tiempo minimo, memoria, asignaciones). Si no lo
# esta, se recurre a una medicion simple basada en `@elapsed` de la
# libreria estandar, de modo que el script funcione sin dependencias.
#
# Para instalar BenchmarkTools (opcional, recomendado):
#     julia -e 'import Pkg; Pkg.add("BenchmarkTools")'

using Printf

include(joinpath(@__DIR__, "..", "src", "numericalMethods.jl"))

# -------------------------------------------------------------------
# Funciones de prueba (mismas raices que en runtests.jl).
# -------------------------------------------------------------------
f1(x)  = x^2 - 2
df1(x) = 2x

f2(x)  = cos(x) - x
df2(x) = -sin(x) - 1

f3(x)  = x^3 - x - 2
df3(x) = 3x^2 - 1

# Detecta si BenchmarkTools esta disponible sin fallar si no lo esta.
const HAS_BENCHMARKTOOLS = try
    @eval using BenchmarkTools
    true
catch
    false
end

"""
    medir(nombre, f)

Mide el tiempo de ejecucion de la closure `f` (sin argumentos) que ejecuta
un metodo. Usa BenchmarkTools si esta disponible; de lo contrario promedia
`@elapsed` sobre varias repeticiones tras un calentamiento (warm-up) del JIT.
"""
function medir(nombre::AbstractString, f)
    if HAS_BENCHMARKTOOLS
        # `$f` interpola la closure para que no se mida el acceso a la global.
        b = @eval @benchmark $f()
        t_ns   = minimum(b.times)   # nanosegundos
        mem    = b.memory           # bytes
        allocs = b.allocs
        @printf("  %-24s %10.3f us  %8d bytes  %5d allocs\n",
                nombre, t_ns / 1e3, mem, allocs)
    else
        f()                         # warm-up (compilacion JIT)
        reps = 10_000
        t = @elapsed for _ in 1:reps
            f()
        end
        @printf("  %-24s %10.3f us  (promedio de %d reps)\n",
                nombre, (t / reps) * 1e6, reps)
    end
end

function main()
    println("="^72)
    println("Benchmark de metodos numericos")
    println(HAS_BENCHMARKTOOLS ? "Backend: BenchmarkTools" :
            "Backend: @elapsed (instala BenchmarkTools para mas precision)")
    println("="^72)

    println("\nf(x) = x^2 - 2   (raiz en sqrt(2))")
    medir("Biseccion [0,2]",   () -> biseccion(f1, 0.0, 2.0))
    medir("Newton (x0=1.0)",   () -> newton(f1, df1, 1.0))
    medir("Secante [0,2]",     () -> secante(f1, 0.0, 2.0))

    println("\nf(x) = cos(x) - x   (raiz en ~0.739)")
    medir("Biseccion [0,1]",   () -> biseccion(f2, 0.0, 1.0))
    medir("Newton (x0=0.5)",   () -> newton(f2, df2, 0.5))
    medir("Secante [0,1]",     () -> secante(f2, 0.0, 1.0))

    println("\nf(x) = x^3 - x - 2   (raiz en ~1.521)")
    medir("Biseccion [1,2]",   () -> biseccion(f3, 1.0, 2.0))
    medir("Newton (x0=1.5)",   () -> newton(f3, df3, 1.5))
    medir("Secante [1,2]",     () -> secante(f3, 1.0, 2.0))

    println("\n" * "="^72)
    println("Iteraciones hasta converger (f(x) = x^2 - 2)")
    println("="^72)
    @printf("  %-12s %d\n", "Biseccion", biseccion(f1, 0.0, 2.0).iteraciones)
    @printf("  %-12s %d\n", "Newton",    newton(f1, df1, 1.0).iteraciones)
    @printf("  %-12s %d\n", "Secante",   secante(f1, 0.0, 2.0).iteraciones)
end

main()
