# src/NumericalMethods.jl


"""
    biseccion(f, a, b; M=100, delta=1e-5, epsilon=1e-5)

Encuentra una raiz de `f` en el intervalo `[a, b]` usando el metodo de Biseccion.
"""
function biseccion(f::Function, a::Float64, b::Float64; 
                   M::Int=100, delta::Float64=1e-5, epsilon::Float64=1e-5)
    u = f(a)
    v = f(b)
    e = b - a
    
    if sign(u) == sign(v)
        error("La funcion debe tener signos opuestos en a y b.")
    end
    
    for k in 1:M
        e = e / 2.0
        c = a + e
        w = f(c)
        
        if abs(e) < delta || abs(w) < epsilon
            return (raiz = c, iteraciones = k, convergio = true)
        end
        
        if sign(w) != sign(u)
            b = c
            v = w
        else
            a = c
            u = w
        end
    end
    
    @warn "Biseccion no convergio en $M iteraciones."
    return (raiz = a + e, iteraciones = M, convergio = false)
end

"""
    newton(f, df, x0; M=50, delta=1e-5, epsilon=1e-5)

Encuentra una raiz de `f` usando su derivada `df` y un punto inicial `x0`.
"""
function newton(f::Function, df::Function, x0::Float64; 
                M::Int=50, delta::Float64=1e-5, epsilon::Float64=1e-5)
    v = f(x0)
    
    if abs(v) < epsilon
        return (raiz = x0, iteraciones = 0, convergio = true)
    end
    
    for k in 1:M
        derivada = df(x0)
        if derivada == 0
            error("Derivada cero en x = $x0")
        end
        
        x1 = x0 - v / derivada
        v = f(x1)
        
        if abs(x1 - x0) < delta || abs(v) < epsilon
            return (raiz = x1, iteraciones = k, convergio = true)
        end
        
        x0 = x1
    end
    
    @warn "Newton no convergio en $M iteraciones."
    return (raiz = x0, iteraciones = M, convergio = false)
end