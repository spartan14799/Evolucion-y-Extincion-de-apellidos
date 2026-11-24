# src/core/graphTopologies.jl
#
# Generadores de topologias para el modelo por parejas. Construyen el
# diccionario de adyacencia `compatibles::Dict{Int,Vector{Int}}` (quien puede
# emparejarse con quien) que consume src/coupleModel.jl. Cada apellido es un
# vertice; por defecto cada apellido es compatible consigo mismo (auto-lazo).
#
# Documentacion completa (topologias disponibles, convenciones, uso):
# docs/code/graphTopologies.md

"""
    _agregar_autolazo!(compatibles)

Asegura que cada vertice sea compatible consigo mismo. Modifica el Dict en sitio.
"""
function _agregar_autolazo!(compatibles::Dict{Int,Vector{Int}})
    for (i, vecinos) in compatibles
        if !(i in vecinos)
            push!(vecinos, i)
        end
        sort!(unique!(vecinos))
    end
    return compatibles
end

"""
    topo_completo(n::Int; autolazo=true) -> Dict{Int,Vector{Int}}

Grafo completo: todos los apellidos 1..n son compatibles con todos. Es la
topologia por defecto del modelo (equivalente a no pasar `compatibles`), pero
explicita. Sirve de baseline.
"""
function topo_completo(n::Int; autolazo::Bool=true)
    apellidos = collect(1:n)
    compatibles = Dict{Int,Vector{Int}}()
    for i in apellidos
        compatibles[i] = autolazo ? copy(apellidos) : filter(!=(i), apellidos)
    end
    return compatibles
end

"""
    topo_cadena(n::Int; autolazo=true) -> Dict{Int,Vector{Int}}

Camino (path) de n apellidos: el apellido `i` es compatible con `i-1`, `i` e
`i+1`. Los extremos (1 y n) tienen un solo vecino ademas de si mismos.

Es la topologia del "apellido puente": en una cadena 1-2-3-4, los apellidos 2 y
3 son puentes entre los extremos. La pregunta de experimento es que pasa si un
puente se extingue antes que los apellidos que conecta.
"""
function topo_cadena(n::Int; autolazo::Bool=true)
    compatibles = Dict{Int,Vector{Int}}()
    for i in 1:n
        vecinos = Int[]
        i > 1 && push!(vecinos, i - 1)
        push!(vecinos, i)
        i < n && push!(vecinos, i + 1)
        compatibles[i] = vecinos
    end
    autolazo || foreach(i -> filter!(!=(i), compatibles[i]), 1:n)
    return compatibles
end

"""
    topo_anillo(n::Int; autolazo=true) -> Dict{Int,Vector{Int}}

Anillo (ciclo): como la cadena pero cerrado, el apellido 1 y el `n` tambien son
vecinos. Util para comparar con la cadena sin el efecto de los extremos.
"""
function topo_anillo(n::Int; autolazo::Bool=true)
    compatibles = topo_cadena(n; autolazo=autolazo)
    if n >= 3
        push!(compatibles[1], n)
        push!(compatibles[n], 1)
    end
    for i in 1:n
        sort!(unique!(compatibles[i]))
    end
    autolazo || foreach(i -> filter!(!=(i), compatibles[i]), 1:n)
    return compatibles
end

"""
    topo_estrella(n::Int; centro=1, autolazo=true) -> Dict{Int,Vector{Int}}

Estrella: un apellido `centro` conectado a todos los demas (hojas), y las hojas
solo conectadas al centro. Modela un "apellido puente" explicito y central: todo
el flujo de apellidos pasa por el. Util para ver si el centro tiende a ganar o a
extinguirse primero.
"""
function topo_estrella(n::Int; centro::Int=1, autolazo::Bool=true)
    @assert 1 <= centro <= n "centro debe estar en 1..n"
    compatibles = Dict{Int,Vector{Int}}()
    hojas = filter(!=(centro), collect(1:n))
    compatibles[centro] = copy(collect(1:n))
    for h in hojas
        compatibles[h] = [centro, h]
    end
    if !autolazo
        foreach(i -> filter!(!=(i), compatibles[i]), 1:n)
    end
    for i in 1:n
        sort!(unique!(compatibles[i]))
    end
    return compatibles
end

"""
    topo_grid(filas::Int, cols::Int; autolazo=true) -> Dict{Int,Vector{Int}}

Grid (malla) cuadrada de `filas` x `cols` apellidos. Cada apellido es compatible
con sus vecinos ortogonales (arriba, abajo, izquierda, derecha) y consigo mismo.
Los apellidos se numeran por filas: el de la fila `r` y columna `c` tiene id
`(r-1)*cols + c`.
"""
function topo_grid(filas::Int, cols::Int; autolazo::Bool=true)
    idx(r, c) = (r - 1) * cols + c
    compatibles = Dict{Int,Vector{Int}}()
    for r in 1:filas, c in 1:cols
        i = idx(r, c)
        vecinos = Int[]
        r > 1     && push!(vecinos, idx(r - 1, c))
        r < filas && push!(vecinos, idx(r + 1, c))
        c > 1     && push!(vecinos, idx(r, c - 1))
        c < cols  && push!(vecinos, idx(r, c + 1))
        autolazo && push!(vecinos, i)
        compatibles[i] = sort(unique(vecinos))
    end
    return compatibles
end

"""
    topo_bipartito(grupoA, grupoB; cruzado=true, autolazo=true) -> Dict{Int,Vector{Int}}

Grafo bipartito entre dos conjuntos de apellidos `grupoA` y `grupoB`.

- `cruzado=true` (por defecto): cada apellido solo es compatible con los del
  OTRO grupo (bipartito clasico). Modela dos familias que solo se reproducen con
  la otra. En este modo el `autolazo` no aplica a la relacion dentro del grupo,
  pero se puede permitir el auto-lazo individual con `autolazo=true`.
- `cruzado=false`: cada grupo es internamente completo y NO hay aristas entre
  grupos. Esto produce dos componentes aisladas (caso Romeo y Julieta): cada
  grupo fija su propio apellido y sobreviven dos linajes de forma permanente.
"""
function topo_bipartito(grupoA::Vector{Int}, grupoB::Vector{Int};
                        cruzado::Bool=true, autolazo::Bool=true)
    compatibles = Dict{Int,Vector{Int}}()
    if cruzado
        for i in grupoA
            compatibles[i] = autolazo ? vcat(grupoB, i) : copy(grupoB)
        end
        for j in grupoB
            compatibles[j] = autolazo ? vcat(grupoA, j) : copy(grupoA)
        end
    else
        for i in grupoA
            compatibles[i] = copy(grupoA)
        end
        for j in grupoB
            compatibles[j] = copy(grupoB)
        end
    end
    for (i, v) in compatibles
        sort!(unique!(v))
    end
    return compatibles
end

"""
    topo_tripartito(gA, gB, gC; ciclico=true, autolazo=true) -> Dict{Int,Vector{Int}}

Grafo tripartito entre tres conjuntos de apellidos.

- `ciclico=true` (por defecto): A se empareja con B, B con C, C con A (anillo de
  grupos). Cada grupo tiene un vecino "aguas arriba" y otro "aguas abajo".
- `ciclico=false`: A-B-C en cadena (A con B, B con A y C, C con B). El grupo B
  actua de puente entre A y C.

`autolazo=true` permite ademas la reproduccion dentro del propio apellido.
"""
function topo_tripartito(gA::Vector{Int}, gB::Vector{Int}, gC::Vector{Int};
                         ciclico::Bool=true, autolazo::Bool=true)
    compatibles = Dict{Int,Vector{Int}}()
    if ciclico
        for i in gA; compatibles[i] = vcat(gB, gC); end  # A limita con B y con C (anillo)
        for i in gB; compatibles[i] = vcat(gA, gC); end
        for i in gC; compatibles[i] = vcat(gA, gB); end
    else
        for i in gA; compatibles[i] = copy(gB); end        # cadena A-B-C
        for i in gB; compatibles[i] = vcat(gA, gC); end    # B es el puente
        for i in gC; compatibles[i] = copy(gB); end
    end
    if autolazo
        for i in vcat(gA, gB, gC)
            push!(compatibles[i], i)
        end
    end
    for (i, v) in compatibles
        sort!(unique!(v))
    end
    return compatibles
end

"""
    topo_componentes(bloques::Vector{Int}...; autolazo=true) -> Dict{Int,Vector{Int}}

Une varios bloques de apellidos en componentes COMPLETAMENTE AISLADAS (sin
aristas entre bloques). Cada bloque es internamente un grafo completo.

Es la unica forma de obtener un equilibrio real con k apellidos en este modelo:
cada componente fija su propio ganador, de modo que sobreviven exactamente
tantos apellidos como componentes (uno por bloque). Generaliza el caso Romeo y
Julieta (`topo_componentes([1,2], [3,4])` deja 2 linajes permanentes).
"""
function topo_componentes(bloques::Vector{Int}...; autolazo::Bool=true)
    compatibles = Dict{Int,Vector{Int}}()
    for bloque in bloques
        for i in bloque
            compatibles[i] = autolazo ? copy(bloque) : filter(!=(i), bloque)
        end
    end
    for (i, v) in compatibles
        sort!(unique!(v))
    end
    return compatibles
end
