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
    Pkg.activate(joinpath(@__DIR__, "..", ".."))
    using Random, StatsBase, Plots, PlutoUI, Printf

    # Motor del modelo por parejas y utilidades del proyecto.
    include(joinpath(@__DIR__, "..", "..", "src", "coupleModel.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "core", "graphTopologies.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "core", "experimentMetrics.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "core", "visualization.jl"))
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

# ╔═╡ c0000001-0000-0000-0000-000000000100
md"""
---
## Area apilada de proporciones por generacion

Misma informacion que el conteo por apellido, pero en proporciones (cada
generacion suma 1). El area total es constante porque el modelo conserva N; lo
que cambia es el reparto entre apellidos a medida que la diversidad baja.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000101
begin
    _hist_area = resultado.historial
    _ngen_area = length(_hist_area)
    _M_area = conteos_por_gen(_hist_area, apellidos_topo)
    _props_area = _M_area ./ sum(_M_area; dims = 2)
    _labels_area = reshape(["Apellido $a" for a in apellidos_topo], 1, :)
    areaplot(0:(_ngen_area - 1), _props_area;
             label = _labels_area,
             xlabel = "Generacion", ylabel = "Proporcion de la poblacion",
             title = "Distribucion de apellidos (area apilada)",
             fillalpha = 0.75, palette = :tab10, legend = :outertopright,
             ylims = (0, 1))
end

# ╔═╡ c0000001-0000-0000-0000-000000000110
md"""
---
## Estado generacion por generacion

Detalle de las primeras y ultimas generaciones de la corrida actual: cuantas
personas tiene cada apellido y cuantos apellidos siguen vivos.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000111
begin
    function _formato_gen_parejas(gen_idx, estado)
        distribucion = sort(collect(contar(estado)))
        resumen = join(["ap $k: $v" for (k, v) in distribucion], ", ")
        n_unicos = length(unique(estado))
        "**Gen $(gen_idx - 1)** ($n_unicos apellidos) -> [$resumen]"
    end

    _hist_det = resultado.historial
    _ng_det = length(_hist_det)
    _primeras = min(5, _ng_det)
    _ult_start = max(_primeras + 1, _ng_det - 4)

    _lineas = String[]
    for i in 1:_primeras
        push!(_lineas, _formato_gen_parejas(i, _hist_det[i]))
    end
    if _ult_start > _primeras + 1
        push!(_lineas, "... ($(_ult_start - _primeras - 1) generaciones intermedias) ...")
    end
    for i in _ult_start:_ng_det
        push!(_lineas, _formato_gen_parejas(i, _hist_det[i]))
    end
    Markdown.parse(join(_lineas, "\n\n"))
end

# ╔═╡ c0000001-0000-0000-0000-000000000120
md"""
---
## Matriz de probabilidad de emparejamiento P(i -> j)

Para la poblacion inicial, cada celda `(i, j)` es la probabilidad de que una
persona del apellido `i` se empareje con una del apellido `j`, segun el modo de
probabilidad elegido y las aristas de la topologia. Las filas suman 1. Las
celdas en 0 son apellidos que **no pueden** emparejarse (sin arista). Es el
analogo, en el modelo por parejas, de la matriz de transicion del grafo.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000121
begin
    _ids_mat = apellidos_topo
    _n_mat = length(_ids_mat)
    _vivos_mat = Set(_ids_mat)
    _mat = zeros(_n_mat, _n_mat)
    for (ii, i) in enumerate(_ids_mat)
        destinos, pesos = pesos_emparejamiento(i, compatibles, conteos_ini, _vivos_mat, Symbol(modo_prob))
        s = sum(pesos)
        if s > 0
            for (d, w) in zip(destinos, pesos)
                _mat[ii, findfirst(==(d), _ids_mat)] = w / s
            end
        end
    end
    _etq_mat = ["ap $a" for a in _ids_mat]
    _cmap = cgrad([:white, "#fef08a", "#fdba74", "#f87171", "#ef4444"], [0.0, 0.05, 0.35, 0.7, 1.0])
    _fsize = _n_mat <= 6 ? 10 : (_n_mat <= 10 ? 8 : 6)
    _pmat = heatmap(_etq_mat, _etq_mat, _mat;
                    yflip = true,
                    title = "P(i -> j) emparejamiento ($(modo_prob))",
                    xlabel = "Apellido pareja (j)", ylabel = "Apellido origen (i)",
                    color = _cmap, clims = (0.0, 1.0),
                    size = (560, 480), aspect_ratio = :equal)
    for ii in 1:_n_mat, jj in 1:_n_mat
        annotate!(_pmat, jj, ii, text(@sprintf("%.2f", _mat[ii, jj]), _fsize, :center, :black))
    end
    _pmat
end

# ╔═╡ c0000001-0000-0000-0000-000000000130
md"""
---
## Distribucion del tiempo hasta el estado absorbente

Cada corrida para al fijar un apellido o al alcanzar un equilibrio permanente
(o al tope de generaciones). Repetimos muchas corridas de la configuracion
actual y vemos cuantas generaciones tarda tipicamente. Solo se cuentan las
corridas que terminaron absorbidas; las transitorias (cortadas por el tope) se
reportan aparte.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000131
begin
    _gen_abs = Int[]
    _n_transit = 0
    for _r in 1:n_reps
        _res = correr_experimento(estado_inicial, compatibles;
                                  modo_prob = Symbol(modo_prob),
                                  max_generaciones = max_gen,
                                  semilla = _r)
        if _res.absorbido
            push!(_gen_abs, _res.generaciones)
        elseif _res.transitorio
            _n_transit += 1
        end
    end

    if isempty(_gen_abs)
        md"Ninguna de las $(n_reps) corridas alcanzo estado absorbente dentro del tope ($(max_gen) generaciones)."
    else
        _p_hist = histogram(_gen_abs;
                            bins = max(10, div(length(_gen_abs), 5)),
                            xlabel = "Generaciones hasta absorcion",
                            ylabel = "Frecuencia",
                            title = "Tiempo hasta absorcion ($(length(_gen_abs))/$(n_reps) corridas)",
                            label = "Corridas absorbidas", color = :steelblue,
                            alpha = 0.75, size = (780, 380))
        vline!(_p_hist, [mean(_gen_abs)];
               label = @sprintf("Media = %.1f", mean(_gen_abs)),
               color = :red, lw = 2, linestyle = :dash)
        _p_hist
    end
end

# ╔═╡ c0000001-0000-0000-0000-000000000132
md"""
### Probabilidad de ganar frente a la frecuencia inicial

En el modelo neutral, la probabilidad de que un apellido sea el unico
superviviente de su componente conexa es proporcional a su frecuencia inicial
**dentro de esa componente**. Comparamos la probabilidad empirica de ganar
(`prob_ganar`, de las repeticiones de arriba) contra la frecuencia inicial de
cada apellido. Para una sola componente conexa, ambas deberian coincidir salvo
ruido muestral.
"""

# ╔═╡ c0000001-0000-0000-0000-000000000133
begin
    _tot_ini = length(estado_inicial)
    _freq_ini = Dict(a => get(conteos_ini, a, 0) / _tot_ini for a in apellidos_topo)
    _filas_pg = ["| Apellido | Personas iniciales | Frecuencia inicial | P(unico superviviente) empirica |",
                 "|:--|:--|:--|:--|"]
    for a in apellidos_topo
        _pg = get(rep.prob_ganar, a, 0.0)
        push!(_filas_pg, "| ap $a | $(get(conteos_ini, a, 0)) | $(round(_freq_ini[a], digits=3)) | $(round(_pg, digits=3)) |")
    end
    Markdown.parse(join(_filas_pg, "\n"))
end

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
# ╟─c0000001-0000-0000-0000-000000000100
# ╠═c0000001-0000-0000-0000-000000000101
# ╟─c0000001-0000-0000-0000-000000000110
# ╟─c0000001-0000-0000-0000-000000000111
# ╟─c0000001-0000-0000-0000-000000000120
# ╠═c0000001-0000-0000-0000-000000000121
# ╟─c0000001-0000-0000-0000-000000000130
# ╠═c0000001-0000-0000-0000-000000000131
# ╟─c0000001-0000-0000-0000-000000000132
# ╟─c0000001-0000-0000-0000-000000000133
# ╟─c0000001-0000-0000-0000-000000000050
