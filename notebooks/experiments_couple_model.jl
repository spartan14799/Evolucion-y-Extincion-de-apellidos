### A Pluto.jl notebook ###
# v1.1.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c0000001-0000-0000-0000-000000000001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using Random, StatsBase, Plots, PlutoUI, Printf

    # Motor del modelo por parejas y utilidades del proyecto.
    include(joinpath(@__DIR__, "..", "src", "modeloParejas.jl"))
    include(joinpath(@__DIR__, "..", "src", "graphTopologies.jl"))
    include(joinpath(@__DIR__, "..", "src", "experimentos.jl"))
    include(joinpath(@__DIR__, "..", "src", "visualizacion.jl"))
    md"Entorno y modulos cargados."
end

# ╔═╡ c0000001-0000-0000-0000-000000000002
md"""
# Modelo por parejas: evolucion y extincion de apellidos

## Que es este modelo

Cada apellido es un vertice de un grafo dirigido. Las aristas (`compatibles`)
indican que apellidos pueden emparejarse. La poblacion se parametriza por el
**numero de parejas** `P`, de modo que el total `N = 2P` siempre es par.

Una generacion transforma los `2P` padres en `2P` hijos:

1. **Emparejamiento sin reemplazo.** Se barajan las personas y, en orden
   aleatorio, cada persona suelta busca pareja entre las personas todavia
   sueltas cuyo apellido sea compatible con el suyo. La pareja concreta se
   sortea segun el modo de probabilidad (uniforme por aristas o ponderado por
   conteo).
2. **Herencia.** Por cada pareja se elige al azar (50/50) a uno de los dos
   miembros. Los **dos hijos** de esa pareja toman ese apellido; el otro no se
   transmite en esa pareja.
3. Cada pareja aporta exactamente 2 hijos, asi que el tamano total se conserva.

Si al emparejar queda gente suelta sin ningun compatible disponible, se empareja
de forma forzada aunque sean incompatibles. Cada evento asi es una **mutacion
definida** y se registra.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000003
md"""
## Por que el tamano de cada grupo se conserva

A diferencia de un modelo donde se reconstruyen `N` cupos globales (que deja
derivar a cada grupo hasta fijar uno solo), aqui la unidad reproductiva es la
**pareja**. Si dos personas del apellido A se emparejan entre si (auto-lazo),
sus dos hijos son A: el grupo A no pierde tamano. La diversidad solo puede bajar
cuando se forma una **pareja mixta** (dos apellidos distintos compatibles), pues
ahi un apellido le gana el relevo al otro por el volado de la herencia.

Consecuencia (criterio exacto): si entre los apellidos vivos **no** existe
ninguna arista entre apellidos distintos, ya no pueden formarse parejas mixtas y
el estado es un **equilibrio permanente**. El caso Romeo y Julieta (dos grupos
sin aristas entre si) mantiene dos apellidos para siempre.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000004
md"""
---
## Laboratorio interactivo

Elige una topologia y sus parametros. El notebook construye el grafo, lo dibuja,
corre una simulacion y muestra como evoluciona.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000010
md"""**Topologia**: $(@bind topologia Select([
    "completo" => "Completo",
    "cadena" => "Cadena (apellido puente)",
    "anillo" => "Anillo",
    "estrella" => "Estrella (puente central)",
    "grid" => "Grid cuadrado",
    "bipartito_cruzado" => "Bipartito cruzado",
    "bipartito_aislado" => "Bipartito aislado (Romeo y Julieta)",
    "tripartito_ciclico" => "Tripartito ciclico",
    "tripartito_cadena" => "Tripartito en cadena",
], default="cadena"))"""

# ╔═╡ c0000001-0000-0000-0000-000000000011
md"""**Numero de apellidos** (para completo / cadena / anillo / estrella): $(@bind n_apellidos Slider(2:1:8, default=4, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000012
md"""**Lado del grid** (grid cuadrado de lado x lado): $(@bind lado_grid Slider(2:1:4, default=2, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000013
md"""**Parejas por apellido** (cada pareja son 2 personas): $(@bind parejas_por_apellido Slider(1:1:10, default=3, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000014
md"""**Modo de probabilidad**: $(@bind modo_prob Select([
    "uniforme_aristas" => "Uniforme por aristas (1 / numero de aristas)",
    "ponderado" => "Ponderado por conteo",
]))"""

# ╔═╡ c0000001-0000-0000-0000-000000000015
md"""**Maximo de generaciones**: $(@bind max_gen Slider(50:50:5000, default=1000, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000016
md"""**Semilla**: $(@bind semilla Slider(1:200, default=7, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000020
begin
    # Construccion de la topologia elegida.
    compatibles = if topologia == "completo"
        topo_completo(n_apellidos)
    elseif topologia == "cadena"
        topo_cadena(n_apellidos)
    elseif topologia == "anillo"
        topo_anillo(n_apellidos)
    elseif topologia == "estrella"
        topo_estrella(n_apellidos; centro=1)
    elseif topologia == "grid"
        topo_grid(lado_grid, lado_grid)
    elseif topologia == "bipartito_cruzado"
        mitad = cld(n_apellidos, 2)
        topo_bipartito(collect(1:mitad), collect(mitad+1:n_apellidos); cruzado=true)
    elseif topologia == "bipartito_aislado"
        mitad = cld(n_apellidos, 2)
        topo_bipartito(collect(1:mitad), collect(mitad+1:n_apellidos); cruzado=false)
    elseif topologia == "tripartito_ciclico"
        topo_tripartito([1], [2], [3])
    elseif topologia == "tripartito_cadena"
        topo_tripartito([1], [2], [3]; ciclico=false)
    else
        topo_completo(n_apellidos)
    end

    apellidos_topo = sort(collect(keys(compatibles)))
    # Estado inicial: mismas parejas por apellido para todos los apellidos de la topologia.
    estado_inicial = poblacion_desde_parejas(Dict(a => parejas_por_apellido for a in apellidos_topo))
    conteos_ini = contar(estado_inicial)

    md"""
    ### Configuracion construida

    - Apellidos: $(length(apellidos_topo)) -> $(apellidos_topo)
    - Poblacion total N = $(length(estado_inicial)) ($(parejas_por_apellido) parejas por apellido)
    - Modo de probabilidad: `$(modo_prob)`
    """
end

# ╔═╡ c0000001-0000-0000-0000-000000000021
md"### Grafo de compatibilidad"

# ╔═╡ c0000001-0000-0000-0000-000000000022
dibujar_grafo(compatibles; conteos=conteos_ini,
              titulo="Topologia: $(topologia)")

# ╔═╡ c0000001-0000-0000-0000-000000000030
begin
    parametros = ParametrosParejas(max_generaciones=max_gen, modo_prob=Symbol(modo_prob))
    resultado = simular_parejas(estado_inicial, compatibles;
                                parametros=parametros, rng=MersenneTwister(semilla))
    resumen = resumen_experimento(resultado)

    estado_label = resultado.absorbido ?
        (resumen.n_vivos_final == 1 ? "Fijacion de un apellido" : "Equilibrio permanente de $(resumen.n_vivos_final) apellidos") :
        (resultado.transitorio ? "Transitorio (corte por tope de generaciones)" : "Detenido")

    md"""
    ### Resultado de la simulacion

    | Campo | Valor |
    |:--|:--|
    | Estado final | **$(estado_label)** |
    | Generaciones | **$(resultado.generaciones)** |
    | Apellidos vivos al final | **$(resumen.vivos_final)** |
    | Mutaciones (parejas forzadas) | **$(resumen.mutaciones)** |
    """
end

# ╔═╡ c0000001-0000-0000-0000-000000000031
md"### Conteo de cada apellido por generacion"

# ╔═╡ c0000001-0000-0000-0000-000000000032
graficar_conteos(resultado.historial; apellidos=apellidos_topo)

# ╔═╡ c0000001-0000-0000-0000-000000000033
md"### Numero de apellidos vivos por generacion"

# ╔═╡ c0000001-0000-0000-0000-000000000034
graficar_diversidad(resultado.historial)

# ╔═╡ c0000001-0000-0000-0000-000000000035
md"### Grafo final (solo apellidos supervivientes)"

# ╔═╡ c0000001-0000-0000-0000-000000000036
begin
    vivos_final = Set(resultado.apellidos_final)
    comp_final = Dict(i => [j for j in compatibles[i] if j in vivos_final]
                      for i in keys(compatibles) if i in vivos_final)
    dibujar_grafo(comp_final; conteos=contar(resultado.estado_final),
                  titulo="Supervivientes: $(sort(collect(vivos_final)))")
end

# ╔═╡ c0000001-0000-0000-0000-000000000040
md"""
---
## Experimentos agregados (muchas corridas)

Una sola corrida es aleatoria. Para sacar conclusiones repetimos el experimento
con muchas semillas y agregamos estadisticas: con que frecuencia se fija un solo
apellido, cuantos apellidos sobreviven, que apellido gana, y (si marcas un
apellido puente) con que frecuencia el puente se extingue primero.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000041
md"""**Numero de repeticiones**: $(@bind n_reps Slider(20:20:500, default=200, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000042
md"""**Apellido puente a vigilar** (0 = ninguno): $(@bind puente_sel Slider(0:1:8, default=0, show_value=true))"""

# ╔═╡ c0000001-0000-0000-0000-000000000043
begin
    puente = puente_sel == 0 ? nothing : puente_sel
    rep = correr_repeticiones(estado_inicial, compatibles;
                              modo_prob=Symbol(modo_prob),
                              max_generaciones=max_gen,
                              n_reps=n_reps,
                              semilla_base=1,
                              puente=puente)

    puente_txt = rep.prob_puente_primero === nothing ? "n/a" :
        @sprintf("%.1f%%", 100 * rep.prob_puente_primero)

    md"""
    ### Resumen de $(n_reps) corridas

    | Metrica | Valor |
    |:--|:--|
    | Prob. de alcanzar estado absorbente | $(@sprintf("%.1f%%", 100*rep.prob_absorcion)) |
    | Prob. de terminar transitorio | $(@sprintf("%.1f%%", 100*rep.prob_transitorio)) |
    | Generaciones (media +/- sd) | $(@sprintf("%.1f +/- %.1f", rep.gen_media, rep.gen_std)) |
    | Mutaciones por corrida (media) | $(@sprintf("%.2f", rep.mutaciones_media)) |
    | Puente ($(puente === nothing ? "ninguno" : puente)) se extingue primero | $(puente_txt) |
    """
end

# ╔═╡ c0000001-0000-0000-0000-000000000044
md"### Distribucion del numero de apellidos supervivientes"

# ╔═╡ c0000001-0000-0000-0000-000000000045
graficar_dist_n_vivos(rep)

# ╔═╡ c0000001-0000-0000-0000-000000000046
md"### Probabilidad de ser el unico superviviente (solo corridas que fijaron uno)"

# ╔═╡ c0000001-0000-0000-0000-000000000047
graficar_prob_ganar(rep)

# ╔═╡ c0000001-0000-0000-0000-000000000050
md"""
---
## Guia de experimentos sugeridos

- **Romeo y Julieta**: topologia `Bipartito aislado`. Dos grupos sin aristas
  entre si. El sistema detecta el equilibrio de entrada y conserva los dos
  apellidos para siempre.
- **Apellido puente (cadena)**: topologia `Cadena` con 4 apellidos. Marca el
  apellido puente = 2 en el selector de puente. Observa con que frecuencia el
  puente se extingue primero y cuantos apellidos quedan al final (a menudo dos
  no adyacentes coexisten).
- **Puente central (estrella)**: topologia `Estrella`, puente = 1. El centro
  esta sobre-expuesto a parejas mixtas pero, cuando sobrevive, tiende a dominar.
- **Grid**: topologia `Grid cuadrado`. Estructura local de vecindad; compara el
  tiempo de fijacion con el grafo completo del mismo numero de apellidos.
- **Modo de probabilidad**: compara `uniforme por aristas` frente a `ponderado`
  para ver como cambia que apellido tiende a ganar.
"""

# ╔═╡ Cell order:
# ╟─c0000001-0000-0000-0000-000000000001
# ╟─c0000001-0000-0000-0000-000000000002
# ╟─c0000001-0000-0000-0000-000000000003
# ╟─c0000001-0000-0000-0000-000000000004
# ╟─c0000001-0000-0000-0000-000000000010
# ╟─c0000001-0000-0000-0000-000000000011
# ╟─c0000001-0000-0000-0000-000000000012
# ╟─c0000001-0000-0000-0000-000000000013
# ╟─c0000001-0000-0000-0000-000000000014
# ╟─c0000001-0000-0000-0000-000000000015
# ╟─c0000001-0000-0000-0000-000000000016
# ╟─c0000001-0000-0000-0000-000000000020
# ╟─c0000001-0000-0000-0000-000000000021
# ╠═c0000001-0000-0000-0000-000000000022
# ╟─c0000001-0000-0000-0000-000000000030
# ╟─c0000001-0000-0000-0000-000000000031
# ╠═c0000001-0000-0000-0000-000000000032
# ╟─c0000001-0000-0000-0000-000000000033
# ╠═c0000001-0000-0000-0000-000000000034
# ╟─c0000001-0000-0000-0000-000000000035
# ╠═c0000001-0000-0000-0000-000000000036
# ╟─c0000001-0000-0000-0000-000000000040
# ╟─c0000001-0000-0000-0000-000000000041
# ╟─c0000001-0000-0000-0000-000000000042
# ╟─c0000001-0000-0000-0000-000000000043
# ╟─c0000001-0000-0000-0000-000000000044
# ╠═c0000001-0000-0000-0000-000000000045
# ╟─c0000001-0000-0000-0000-000000000046
# ╠═c0000001-0000-0000-0000-000000000047
# ╟─c0000001-0000-0000-0000-000000000050
