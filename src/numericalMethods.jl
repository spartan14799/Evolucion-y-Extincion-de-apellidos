"""
    biseccion(f, a, b; M=100, delta=1e-5, epsilon=1e-5)

Encuentra una raíz de `f` en el intervalo `[a, b]` usando el método de Bisección.
"""
function biseccion(f::Function, a::Float64, b::Float64; 
                   M::Int=100, delta::Float64=1e-5, epsilon::Float64=1e-5)
    u = f(a)
    v = f(b)
    e = b - a
    
    if sign(u) == sign(v)
        error("La función debe tener signos opuestos en a y b.")
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
    
    @warn "Bisección no convergió en $M iteraciones."
    return (raiz = a + e, iteraciones = M, convergio = false)
end