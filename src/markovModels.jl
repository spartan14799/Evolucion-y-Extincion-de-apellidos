using Distributions
using Random

"""
    ModeloPBGW

Estructura que representa un Proceso de Bienayme-Galton-Watson para unico apellidos.
Almacena la distribucion de probabilidad de descendencia.
"""
struct ModeloPBGW
    pk::Vector{Float64}

    function ModeloPBGW(pk::Vector{Float64})
        if any(p -> p < 0.0, pk) error("Probabilidades negativas.") end
        if !isapprox(sum(pk), 1.0, atol=1e-5) error("No suman 1.0") end
        new(pk)
    end

    function ModeloPBGW(tasa_natalidad::Float64; max_hijos::Int=15)
        # Poisson lambda = tasa_natalidad
        dist = Poisson(tasa_natalidad)
        
        pk = [pdf(dist, k) for k in 0:max_hijos]
        
        # Normalizamos para segurar que sumen exactamente 1
        pk = pk ./ sum(pk)
    
        new(pk)
    end
end

"""
Calcula el promedio de hijos por individuo (Esperanza m = E[X]).
"""
function esperanza(modelo::ModeloPBGW)
    return sum((k - 1) * modelo.pk[k] for k in 1:length(modelo.pk))
end

"""
Verifica si la extincion es matematicamente 100% segura (m <= 1).
"""
function es_extincion_segura(modelo::ModeloPBGW)
    return esperanza(modelo) <= 1.0
end

"""
Evalua la Funcion Generadora de Probabilidad
"""
function fgp(modelo::ModeloPBGW, t::Float64)
    return sum(modelo.pk[k] * (t^(k - 1)) for k in 1:length(modelo.pk))
end



"""
    probabilidad_extincion(modelo; buscador=:biseccion)

Devuelve un `Float64` en [0, 1].
"""
function probabilidad_extincion(modelo::ModeloPBGW; buscador::Symbol=:biseccion)
    # Caso subcritico o critico: la extincion es segura.
    if es_extincion_segura(modelo)
        return 1.0
    end

    # Caso supercritico (m > 1): buscamos el punto fijo q en (0, 1).
    # h(t) = g(t) - t.  Se cumple h(0) = p0 > 0  y  h(1) = 0.
    # El menor punto fijo esta en (0, 1); buscamos en [0, 1 - margen]
    # para evitar capturar la raiz trivial t = 1.
    h(t::Float64) = fgp(modelo, t) - t

    margen = 1e-6
    a, b = 0.0, 1.0 - margen

    # Si en [a, b] no hay cambio de signo, el punto fijo esta muy cerca de 1.
    if sign(h(a)) == sign(h(b))
        return 1.0
    end

    resultado = if buscador == :biseccion
        biseccion(h, a, b)
    elseif buscador == :secante
        # La secante es un metodo ABIERTO: no queda confinado al intervalo.
        # Como t = 1 siempre es raiz de h (la raiz trivial), arrancar cerca de
        # 1 hace que la secante converja a esa raiz en vez de a q en (0, 1).
        # Para evitarlo, primero acotamos q con unos pocos pasos de biseccion
        # y luego lanzamos la secante desde ese intervalo estrecho, donde el
        # unico punto fijo es q.
        pre = biseccion(h, a, b; M = 20)
        c = pre.raiz
        secante(h, max(c - 0.05, 0.0), min(c + 0.05, b))
    else
        error("Buscador no soportado: $buscador (usa :biseccion o :secante)")
    end

    q = clamp(resultado.raiz, 0.0, 1.0)
    return q
end

"""
    simular_generaciones(modelo, generaciones; poblacion_inicial=1, rng=Random.default_rng())

Simula la evolucion estocastica del numero de portadores de un apellido a lo
largo de `generaciones` generaciones, partiendo de `poblacion_inicial`.

En cada generacion, cada portador vivo genera un numero aleatorio de hijos
muestreado de la distribucion de descendencia `pk` del modelo. La poblacion
de la siguiente generacion es la suma de todos los hijos.

Devuelve un `Vector{Int}` de longitud `generaciones + 1` con la poblacion en
cada paso (indice 1 = generacion 0 = poblacion inicial). Si la poblacion llega
a 0, el apellido se extingue y permanece en 0.

`tope_poblacion` acota el tamano maximo del linaje. En modelos supercriticos
(m > 1) la poblacion de los linajes que no se extinguen crece de forma
exponencial, lo que haria la simulacion inviable. Al superar el tope, el
linaje se considera "establecido" (ya no se extinguira) y la poblacion se
fija al tope, evitando iterar sobre millones de individuos sin alterar el
resultado de extincion.
"""
function simular_generaciones(modelo::ModeloPBGW, generaciones::Int;
                              poblacion_inicial::Int=1,
                              tope_poblacion::Int=1_000_000,
                              rng::AbstractRNG=Random.default_rng())
    generaciones >= 0 || error("generaciones debe ser >= 0")
    poblacion_inicial >= 0 || error("poblacion_inicial debe ser >= 0")
    tope_poblacion >= 1 || error("tope_poblacion debe ser >= 1")

    # Distribucion categorica: pk[k] es P(k-1 hijos).
    dist_hijos = Categorical(modelo.pk)

    trayectoria = Vector{Int}(undef, generaciones + 1)
    trayectoria[1] = poblacion_inicial
    poblacion = poblacion_inicial

    for g in 1:generaciones
        if poblacion == 0
            trayectoria[g + 1] = 0
            continue
        end
        # Si el linaje ya supero el tope, se considera establecido: no se
        # extinguira, asi que mantenemos la poblacion en el tope.
        if poblacion >= tope_poblacion
            trayectoria[g + 1] = tope_poblacion
            poblacion = tope_poblacion
            continue
        end
        # Cada individuo produce (muestra - 1) hijos, porque la categoria 1
        # corresponde a 0 hijos.
        hijos_totales = 0
        for _ in 1:poblacion
            hijos_totales += rand(rng, dist_hijos) - 1
        end
        poblacion = hijos_totales
        trayectoria[g + 1] = poblacion
    end

    return trayectoria
end

"""
    frecuencia_extincion(modelo, generaciones; n_sim=10_000, poblacion_inicial=1, rng=...)

Estima por Monte Carlo la fraccion de linajes que se extinguen dentro de
`generaciones` generaciones, ejecutando `n_sim` simulaciones independientes.

Devuelve un `Float64` en [0, 1]: proporcion de simulaciones que terminaron
con poblacion 0. Util para contrastar contra `probabilidad_extincion`.
"""
function frecuencia_extincion(modelo::ModeloPBGW, generaciones::Int;
                              n_sim::Int=10_000,
                              poblacion_inicial::Int=1,
                              tope_poblacion::Int=1_000_000,
                              rng::AbstractRNG=Random.default_rng())
    extintos = 0
    for _ in 1:n_sim
        traj = simular_generaciones(modelo, generaciones;
                                    poblacion_inicial=poblacion_inicial,
                                    tope_poblacion=tope_poblacion, rng=rng)
        if traj[end] == 0
            extintos += 1
        end
    end
    return extintos / n_sim
end
