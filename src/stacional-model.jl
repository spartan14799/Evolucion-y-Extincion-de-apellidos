# src/stacional-model.jl
#
# Modelo "estacionario" (poblacion de tamano constante) para la evolucion y
# extincion de apellidos.
#
# Idea del modelo
# ---------------
# La poblacion se mantiene constante en cada generacion. En cada paso los
# individuos se emparejan al azar (parejas de tamano fijo = 2). De cada pareja
# sobrevive/gana UN apellido, que se propaga a los dos descendientes de esa
# pareja. Asi el tamano total no cambia, pero la diversidad de apellidos solo
# puede mantenerse o disminuir: el proceso avanza siempre hacia la FIJACION de
# un unico apellido (estado absorbente), analogo a la deriva genetica.
#
# Parametros configurables de la dinamica estan reunidos en `ParametrosModelo`.

using Random

"""
    ParametrosModelo

Agrupa los parametros configurables de la dinamica estacionaria.

Campos
- `tamano_pareja::Int`: numero de individuos que se agrupan para elegir un
  unico apellido ganador (por defecto `2`). El tamano de la poblacion debe ser
  divisible por este valor para conservarse constante.
- `max_generaciones::Int`: tope de seguridad de generaciones simuladas antes de
  rendirse si no se alcanza la fijacion (por defecto `10_000`).

Ejemplo
```julia
p = ParametrosModelo(tamano_pareja = 2, max_generaciones = 5_000)
```
"""
Base.@kwdef struct ParametrosModelo
    tamano_pareja::Int = 2
    max_generaciones::Int = 10_000
end

"""
    siguiente_generacion(poblacion, [rng]; parametros = ParametrosModelo())

Simula una transicion de estado en la poblacion bajo el modelo de parejas de
tamano constante.

Mezcla la poblacion al azar, la agrupa en bloques de `parametros.tamano_pareja`
individuos y propaga el apellido del primer miembro de cada bloque a todos los
miembros de ese bloque. Conserva el tamano total de la poblacion.

Lanza un error si `length(poblacion)` no es divisible por `tamano_pareja`, ya
que en ese caso el tamano de la poblacion no se conservaria.
"""
function siguiente_generacion(poblacion::Vector{Int},
                              rng::AbstractRNG = Random.default_rng();
                              parametros::ParametrosModelo = ParametrosModelo())
    k = parametros.tamano_pareja
    n = length(poblacion)
    if k < 1
        error("El tamano de pareja debe ser >= 1 (recibido $k).")
    end
    if n % k != 0
        error("La poblacion ($n) debe ser divisible por el tamano de pareja ($k) para conservarse constante.")
    end

    mezcla = shuffle(rng, poblacion)
    nueva_poblacion = Vector{Int}(undef, n)

    for i in 1:k:n
        apellido_ganador = mezcla[i]
        for j in 0:(k - 1)
            nueva_poblacion[i + j] = apellido_ganador
        end
    end

    return nueva_poblacion
end

"""
    simular_fijacion(estado_inicial; max_generaciones, rng, parametros)

Simula la evolucion hasta la fijacion de un unico apellido, o hasta alcanzar el
limite de generaciones.

Argumentos
- `estado_inicial::Vector{Int}`: apellidos (ids enteros) de la poblacion inicial.

Palabras clave
- `parametros::ParametrosModelo`: parametros de la dinamica. Por compatibilidad,
  `max_generaciones` puede pasarse directamente y sobreescribe el valor de
  `parametros`.
- `max_generaciones::Int`: tope de generaciones (por defecto el de `parametros`).
- `rng::AbstractRNG`: generador de numeros aleatorios (para reproducibilidad).

Devuelve una named tuple con:
- `ganador`: id del apellido fijado, o `nothing` si no hubo fijacion.
- `generaciones`: numero de transiciones simuladas.
- `historial`: vector con el estado de la poblacion en cada generacion.
- `fijacion`: `true` si se alcanzo un unico apellido.
"""
function simular_fijacion(estado_inicial::Vector{Int};
                          parametros::ParametrosModelo = ParametrosModelo(),
                          max_generaciones::Int = parametros.max_generaciones,
                          rng::AbstractRNG = Random.default_rng())
    poblacion = copy(estado_inicial)
    historial = [copy(poblacion)]
    generaciones = 0

    # Continua mientras haya mas de 1 apellido Y no choquemos con el limite.
    while length(unique(poblacion)) > 1 && generaciones < max_generaciones
        poblacion = siguiente_generacion(poblacion, rng; parametros = parametros)
        push!(historial, copy(poblacion))
        generaciones += 1
    end

    # Determinamos si el ciclo termino porque gano alguien, o porque se acabo
    # el tiempo.
    logro_fijacion = length(unique(poblacion)) == 1
    ganador = logro_fijacion ? poblacion[1] : nothing

    return (ganador = ganador, generaciones = generaciones, historial = historial, fijacion = logro_fijacion)
end
