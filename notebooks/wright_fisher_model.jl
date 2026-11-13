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

# ╔═╡ 00000001-0001-0001-0001-000000000001
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using Random, StatsBase, Plots, PlutoUI, Printf

	# ── Cargamos los módulos del proyecto tal cual están en src/ ──
	include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))
	# Esto nos da: GrafoApellidos, siguiente_generacion_grafo, ParametrosModelo, simular_fijacion
	# (stacional-model.jl incluye internamente graphModel.jl)
end

# ╔═╡ 00000001-0001-0001-0001-000000000002
md"""
# 🏛️ Modelo Wright-Fisher para Apellidos
## Evolución y Fijación con Población Constante

---

### 📖 ¿Qué es este modelo?

El **modelo de Wright-Fisher** (originalmente diseñado para genética de poblaciones) modela cómo evoluciona la frecuencia de variantes en una **población de tamaño fijo**. En nuestro caso, los "alelos" son **apellidos**.

La idea central es (tal como la implementa [`src/stacional-model.jl`](../src/stacional-model.jl) con [`src/graphModel.jl`](../src/graphModel.jl)):

> En cada generación, se construye un **grafo dirigido de compatibilidad** entre apellidos (`GrafoApellidos`). Para producir los `N` hijos de la siguiente generación:
> 1. Se elige un apellido **origen** `i` ponderado por su frecuencia en la población actual.
> 2. Se elige una **pareja** `j` según la distribución de probabilidad saliente `P(i → j)` del grafo.
> 3. El hijo hereda el apellido de `i` con probabilidad `p_herencia`, o el de `j` con probabilidad `1 - p_herencia`.

Esto es lo que hace la función `siguiente_generacion_grafo()` de `graphModel.jl`.

El proceso avanza hacia la **fijación** de un solo apellido (estado absorbente), análogo a la **deriva genética**.

!!! info "Analogía biológica"
    Es el mismo fenómeno por el cual, sin selección natural, una mutación neutra termina fijándose o perdiéndose en una población finita. Aquí los "genes" son apellidos.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000003
md"""
### ⛓️ La Cadena de Markov

El estado del sistema en cada generación es el **vector de apellidos** de todos los individuos. Como la transición sólo depende del estado actual (muestreo ponderado + herencia estocástica), el proceso es una **cadena de Markov**.

**Propiedades clave de `simular_fijacion()` (stacional-model.jl):**
- **Espacio de estados**: Todas las posibles asignaciones de apellidos a `N` individuos
- **Transiciones**: `siguiente_generacion_grafo()` → construir grafo → muestrear origen → elegir pareja → heredar apellido
- **Estado absorbente**: Un único apellido en toda la población (`length(unique(poblacion)) == 1`)
- **Irreversibilidad**: La diversidad **nunca aumenta** → la fijación es **inevitable**
- **Invariante**: El tamaño de la población se conserva en cada paso (`N` constante)

La pregunta no es **si** un apellido domina, sino **cuándo** y **cuál**.

**Novedad — Modelo de Grafos:**
El grafo de compatibilidad (`GrafoApellidos`) permite definir qué apellidos pueden interactuar entre sí. En un **grafo completo** (por defecto), todos compiten con todos. Pero podemos definir redes **restringidas** (tipo "mundo pequeño") donde cada apellido solo interactúa con ciertos vecinos, ¡cambiando drásticamente la dinámica de fijación!
"""

# ╔═╡ 00000001-0001-0001-0001-000000000004
md"""
### 🎛️ Panel de Control

Ajusta los parámetros de la simulación:
"""

# ╔═╡ 00000001-0001-0001-0001-000000000005
md"**Población total** (N): $(@bind N_slider Slider(4:2:60, default=10, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000006
md"**Número de apellidos iniciales**: $(@bind num_apellidos_slider Slider(2:1:20, default=4, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000007
md"**Semilla aleatoria** (para reproducibilidad): $(@bind semilla_slider Slider(1:100, default=42, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000008
md"**Máximo de generaciones**: $(@bind max_gen_slider Slider(100:100:5000, default=1000, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000040
md"**Probabilidad de herencia** (`p_herencia` — prob. de heredar el apellido del origen): $(@bind p_herencia_slider Slider(0.0:0.05:1.0, default=0.5, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000009
begin
	# Ajustar apellidos al máximo posible
	num_apellidos = min(num_apellidos_slider, N_slider)
	
	rng_sim = MersenneTwister(semilla_slider)

	# ── Construimos la población inicial ──
	# Aseguramos que cada apellido (1:num_apellidos) aparezca al menos una vez
	# El resto se rellena al azar y se mezcla
	estado_inicial = begin
		pob = collect(1:num_apellidos)
		faltantes = N_slider - num_apellidos
		if faltantes > 0
			append!(pob, rand(rng_sim, 1:num_apellidos, faltantes))
		end
		shuffle(rng_sim, pob)
	end
	
	# ── Usamos ParametrosModelo y simular_fijacion del proyecto ──
	parametros = ParametrosModelo(
		max_generaciones = max_gen_slider,
		p_herencia = p_herencia_slider
	)
	resultado = simular_fijacion(estado_inicial; parametros=parametros, rng=rng_sim)
	
	md"""
	### 📊 Resultado de la Simulación
	
	| Parámetro | Valor |
	|:--|:--|
	| Población (N) | **$(N_slider)** |
	| Apellidos iniciales | **$(num_apellidos)** |
	| p\_herencia | **$(p_herencia_slider)** |
	| Semilla | **$(semilla_slider)** |
	| ¿Hubo fijación? | **$(resultado.fijacion ? "✅ SÍ" : "❌ NO (límite alcanzado)")** |
	| Generaciones transcurridas | **$(resultado.generaciones)** |
	| Apellido ganador | **$(resultado.fijacion ? "ID $(resultado.ganador)" : "Ninguno aún")** |
	"""
end

# ╔═╡ 00000001-0001-0001-0001-000000000010
md"""
### 📈 Evolución de la Frecuencia de Apellidos por Generación

Esta gráfica muestra cómo las proporciones de cada apellido cambian generación tras generación, hasta que uno de ellos alcanza el 100% (fijación).

Los datos provienen de `resultado.historial`, que es el vector de estados devuelto por `simular_fijacion()`.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000011
begin
	# Usamos resultado.historial directamente (output de simular_fijacion)
	historial = resultado.historial
	n_gens = length(historial)
	apellidos_unicos = sort(unique(estado_inicial))
	n_ap = length(apellidos_unicos)
	
	# Calcular frecuencias relativas a partir del historial
	frecuencias = zeros(n_gens, n_ap)
	for (g, estado) in enumerate(historial)
		conteo = countmap(estado)
		for (j, ap) in enumerate(apellidos_unicos)
			frecuencias[g, j] = get(conteo, ap, 0) / length(estado)
		end
	end
	
	colores = palette(:tab10, n_ap)
	
	p1 = plot(
		size = (800, 450),
		xlabel = "Generación",
		ylabel = "Frecuencia relativa",
		title = "Evolución de Frecuencias de Apellidos",
		legend = :outerright,
		ylims = (0, 1.05),
		grid = true, gridstyle = :dash, gridalpha = 0.3
	)
	
	for (j, ap) in enumerate(apellidos_unicos)
		plot!(p1, 0:(n_gens-1), frecuencias[:, j],
			label = "Apellido $ap", linewidth = 2.5,
			color = colores[j], alpha = 0.85)
	end
	
	if resultado.fijacion
		vline!(p1, [resultado.generaciones], 
			label = "Fijación (gen $(resultado.generaciones))",
			color = :red, linestyle = :dash, linewidth = 1.5)
	end
	
	p1
end

# ╔═╡ 00000001-0001-0001-0001-000000000012
md"""
### 📊 Gráfica de Áreas Apiladas

El área total siempre es 1 (100% de la población se conserva — invariante de `siguiente_generacion_grafo()`).
"""

# ╔═╡ 00000001-0001-0001-0001-000000000013
begin
	p2 = areaplot(
		0:(n_gens-1), frecuencias,
		label = reshape(["Apellido $a" for a in apellidos_unicos], 1, :),
		xlabel = "Generación", ylabel = "Proporción de la población",
		title = "Distribución de Apellidos (Áreas Apiladas)",
		size = (800, 400), fillalpha = 0.7, palette = :tab10,
		legend = :outerright
	)
	p2
end

# ╔═╡ 00000001-0001-0001-0001-000000000014
md"""
### 🔬 Estado Generación por Generación

Detalle de las primeras y últimas generaciones:
"""

# ╔═╡ 00000001-0001-0001-0001-000000000015
begin
	# Formato similar al que usa el script: "Gen X | Frecuencias: [ID k:v, ...]"
	function formato_gen(gen_idx, estado)
		n_gen = gen_idx - 1
		distribucion = countmap(estado)
		resumen = join(["ID $k:$v" for (k, v) in sort(collect(distribucion))], ", ")
		n_unicos = length(unique(estado))
		"**Gen $n_gen** ($n_unicos apellidos) → [$resumen]"
	end
	
	primeras = min(5, n_gens)
	ultimas_start = max(1, n_gens - 4)
	
	lineas = String[]
	for i in 1:primeras
		push!(lineas, formato_gen(i, historial[i]))
	end
	if n_gens > 10
		push!(lineas, "⋮ *($(n_gens - 10) generaciones intermedias)* ⋮")
	end
	if ultimas_start > primeras
		for i in ultimas_start:n_gens
			push!(lineas, formato_gen(i, historial[i]))
		end
	end
	
	Markdown.parse(join(lineas, "\n\n"))
end

# ╔═╡ 00000001-0001-0001-0001-000000000041
md"""
---

## 🕸️ Visualización del Grafo de Compatibilidad

El modelo de grafos (`GrafoApellidos`) define **qué apellidos pueden interactuar entre sí**. Por defecto es un **grafo completo** (todos con todos), pero podemos restringirlo para simular redes de **mundo pequeño**.

Veamos cómo se ve el grafo y la matriz de probabilidades de transición para la población actual:
"""

# ╔═╡ 00000001-0001-0001-0001-000000000042
md"**Tipo de grafo:** $(@bind tipo_grafo Select([\"completo\" => \"🔗 Grafo Completo (todos interactúan con todos)\", \"anillo\" => \"⭕ Anillo (cada apellido solo con vecinos ±1)\", \"mundo_pequeno\" => \"🌐 Mundo Pequeño (anillo + atajos aleatorios)\"], default=\"completo\"))"

# ╔═╡ 00000001-0001-0001-0001-000000000043
begin
	# Construir la compatibilidad según el tipo de grafo seleccionado
	apellidos_en_pob = sort(unique(estado_inicial))
	n_apellidos_pob = length(apellidos_en_pob)
	
	function construir_compatibles(aps, tipo, rng_grafo)
		comp = Dict{Int,Vector{Int}}()
		n = length(aps)
		if tipo == "completo"
			for a in aps
				comp[a] = copy(aps)
			end
		elseif tipo == "anillo"
			for (idx, a) in enumerate(aps)
				vecino_izq = aps[mod1(idx - 1, n)]
				vecino_der = aps[mod1(idx + 1, n)]
				comp[a] = sort(unique([a, vecino_izq, vecino_der]))
			end
		elseif tipo == "mundo_pequeno"
			# Anillo base + atajos aleatorios (prob 0.3 de añadir cada no-vecino)
			for (idx, a) in enumerate(aps)
				vecino_izq = aps[mod1(idx - 1, n)]
				vecino_der = aps[mod1(idx + 1, n)]
				vecinos = Set([a, vecino_izq, vecino_der])
				for other in aps
					if other ∉ vecinos && rand(rng_grafo) < 0.3
						push!(vecinos, other)
					end
				end
				comp[a] = sort(collect(vecinos))
			end
		end
		return comp
	end
	
	rng_grafo_vis = MersenneTwister(semilla_slider + 999)
	compatibles_custom = construir_compatibles(apellidos_en_pob, tipo_grafo, rng_grafo_vis)
	
	# Construir el grafo para visualización
	grafo_vis = GrafoApellidos(estado_inicial; compatibles=compatibles_custom)
	
	# ── Visualizar la MATRIZ de probabilidades ──
	n_a = length(grafo_vis.apellidos)
	mat_probs = zeros(n_a, n_a)
	for (i_idx, i) in enumerate(grafo_vis.apellidos)
		probs_i = get(grafo_vis.probabilidades, i, Dict{Int,Float64}())
		for (j_idx, j) in enumerate(grafo_vis.apellidos)
			mat_probs[i_idx, j_idx] = get(probs_i, j, 0.0)
		end
	end
	
	etiquetas_aps = ["Ap. $a" for a in grafo_vis.apellidos]
	
	# Paleta cálida hacia tonos rojos brillantes para valores altos (sin oscurecer el fondo)
	cmap_calor = cgrad([:white, "#fef08a", "#fdba74", "#f87171", "#ef4444"], [0.0, 0.05, 0.35, 0.7, 1.0])
	fsize = n_a <= 6 ? 12 : (n_a <= 10 ? 10 : 8)
	font_fam = Sys.iswindows() ? "arialbd" : "sans-serif"
	
	p_mat = heatmap(
		etiquetas_aps, etiquetas_aps, mat_probs,
		yflip = true,
		title = "Matriz de Transición P(i → j) — Grafo $(tipo_grafo)",
		xlabel = "Apellido destino (j)", ylabel = "Apellido origen (i)",
		color = cmap_calor, clims = (0.0, 1.0),
		size = (620, 520), aspect_ratio = :equal
	)
	
	# Anotar valores en CADA cuadro correspondiente, perfectamente centrados y en NEGRO NEGRILLA
	for i_idx in 1:n_a
		for j_idx in 1:n_a
			val = mat_probs[i_idx, j_idx]
			txt = @sprintf("%.2f", val)
			annotate!(p_mat, j_idx, i_idx, text(txt, font(fsize, :center, :black, font_fam)))
		end
	end
	
	p_mat
end

# ╔═╡ 00000001-0001-0001-0001-000000000044
begin
	# ── Visualizar el grafo como red de nodos ──
	n_nodos = length(grafo_vis.apellidos)
	angulos = range(0, 2π, length=n_nodos + 1)[1:end-1]
	radio = 3.0
	xs_nodos = radio .* cos.(angulos)
	ys_nodos = radio .* sin.(angulos)
	
	p_red = plot(
		size = (600, 600),
		title = "Red de Compatibilidad — Grafo $(tipo_grafo)",
		legend = false,
		aspect_ratio = :equal,
		xlims = (-5, 5), ylims = (-5, 5),
		grid = false, framestyle = :none
	)
	
	# Dibujar aristas (conexiones de compatibilidad)
	for (i_idx, i) in enumerate(grafo_vis.apellidos)
		for j in get(grafo_vis.compatibles, i, Int[])
			j_idx = findfirst(==(j), grafo_vis.apellidos)
			if j_idx !== nothing && j != i
				prob_ij = get(get(grafo_vis.probabilidades, i, Dict{Int,Float64}()), j, 0.0)
				plot!(p_red, 
					[xs_nodos[i_idx], xs_nodos[j_idx]],
					[ys_nodos[i_idx], ys_nodos[j_idx]],
					color = :gray, alpha = 0.3 + 0.7 * prob_ij,
					linewidth = 1 + 3 * prob_ij)
			end
		end
	end
	
	# Dibujar nodos
	node_colors = palette(:tab10, n_nodos)
	conteo_pob = countmap(estado_inicial)
	for (i_idx, ap) in enumerate(grafo_vis.apellidos)
		cnt = get(conteo_pob, ap, 0)
		scatter!(p_red, [xs_nodos[i_idx]], [ys_nodos[i_idx]],
			color = node_colors[i_idx],
			markersize = 12 + 3 * cnt,
			markerstrokecolor = :black,
			markerstrokewidth = 2)
		annotate!(p_red, xs_nodos[i_idx], ys_nodos[i_idx] + 0.7,
			text("Ap. $ap\n(n=$cnt)", 8, :center))
	end
	
	p_red
end

# ╔═╡ 00000001-0001-0001-0001-000000000045
md"""
> **Lectura de las visualizaciones:**
> - El **heatmap** muestra la probabilidad de que un individuo con apellido $i$ (fila) elija como pareja a alguien con apellido $j$ (columna). En un grafo completo, las probabilidades son proporcionales a la frecuencia poblacional. En grafos restringidos, las celdas con 0 indican que esos apellidos **no pueden interactuar**.
> - La **red** muestra los nodos (apellidos) y aristas (compatibilidades). El grosor de las aristas refleja la probabilidad de transición. El tamaño de cada nodo refleja su población.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000046
md"""
---

### 🔄 Simulación con Grafo Personalizado

Ahora corremos `simular_fijacion()` usando el grafo seleccionado arriba para ver cómo la topología de la red afecta la dinámica de fijación:
"""

# ╔═╡ 00000001-0001-0001-0001-000000000047
begin
	rng_grafo_sim = MersenneTwister(semilla_slider)
	
	# Reconstruir población inicial con la misma semilla
	estado_grafo = begin
		pob_g = collect(1:num_apellidos)
		falt_g = N_slider - num_apellidos
		if falt_g > 0
			append!(pob_g, rand(rng_grafo_sim, 1:num_apellidos, falt_g))
		end
		shuffle(rng_grafo_sim, pob_g)
	end
	
	# Construir compatibles para el grafo
	rng_gc = MersenneTwister(semilla_slider + 999)
	compat_sim = construir_compatibles(sort(unique(estado_grafo)), tipo_grafo, rng_gc)
	
	# Función para calcular probabilidades custom a partir de compatibles
	function prob_custom_fn(poblacion)
		g = GrafoApellidos(poblacion; compatibles=compat_sim)
		return g.probabilidades
	end
	
	# Configurar parámetros con probabilidades custom si no es grafo completo
	params_grafo = if tipo_grafo == "completo"
		ParametrosModelo(max_generaciones=max_gen_slider, p_herencia=p_herencia_slider)
	else
		ParametrosModelo(
			max_generaciones=max_gen_slider, 
			p_herencia=p_herencia_slider,
			fijar_probabilidades=prob_custom_fn
		)
	end
	
	rng_sim_grafo = MersenneTwister(semilla_slider + 500)
	resultado_grafo = simular_fijacion(estado_grafo; parametros=params_grafo, rng=rng_sim_grafo)
	
	# Calcular frecuencias
	hist_g = resultado_grafo.historial
	n_gens_g = length(hist_g)
	aps_g = sort(unique(estado_grafo))
	n_ap_g = length(aps_g)
	
	freqs_g = zeros(n_gens_g, n_ap_g)
	for (gen, est) in enumerate(hist_g)
		cnt = countmap(est)
		for (j, ap) in enumerate(aps_g)
			freqs_g[gen, j] = get(cnt, ap, 0) / length(est)
		end
	end
	
	p_grafo_freq = plot(
		size = (800, 450),
		xlabel = "Generación", ylabel = "Frecuencia relativa",
		title = "Evolución con Grafo $(tipo_grafo) (p_herencia=$(p_herencia_slider))",
		legend = :outerright,
		ylims = (0, 1.05),
		grid = true, gridstyle = :dash, gridalpha = 0.3
	)
	
	colores_g = palette(:tab10, n_ap_g)
	for (j, ap) in enumerate(aps_g)
		plot!(p_grafo_freq, 0:(n_gens_g-1), freqs_g[:, j],
			label = "Apellido $ap", linewidth = 2.5,
			color = colores_g[j], alpha = 0.85)
	end
	
	if resultado_grafo.fijacion
		vline!(p_grafo_freq, [resultado_grafo.generaciones],
			label = "Fijación (gen $(resultado_grafo.generaciones))",
			color = :red, linestyle = :dash, linewidth = 1.5)
	end
	
	p_grafo_freq
end

# ╔═╡ 00000001-0001-0001-0001-000000000048
begin
	md"""
	**Resultado con grafo $(tipo_grafo):**
	
	| Parámetro | Valor |
	|:--|:--|
	| ¿Hubo fijación? | **$(resultado_grafo.fijacion ? "✅ SÍ" : "❌ NO")** |
	| Generaciones | **$(resultado_grafo.generaciones)** |
	| Ganador | **$(resultado_grafo.fijacion ? "Apellido $(resultado_grafo.ganador)" : "Ninguno aún")** |
	| Topología | **$(tipo_grafo)** |
	| p\_herencia | **$(p_herencia_slider)** |
	
	> [!TIP]
	> Prueba cambiar entre **grafo completo**, **anillo** y **mundo pequeño** para ver cómo la estructura de la red afecta el tiempo y resultado de la fijación. En grafos restringidos, los apellidos que no están conectados **no pueden competir directamente**, lo que puede retrasar o alterar la fijación.
	"""
end

# ╔═╡ 00000001-0001-0001-0001-000000000030
md"""
---

## 🎮 Minijuego: El Torneo de Linajes de Colombia
### Árbol Genealógico de Coalescencia (Censo DANE 2018)

Aquí simulamos el árbol genealógico completo generación por generación, usando el motor `siguiente_generacion_grafo()` del proyecto:
- **Cada fila horizontal** representa una generación sucesiva.
- **Cada círculo** representa un individuo de la población ($N = 10$). El color indica su apellido.
- **Las líneas negras** muestran la herencia directa: conectan a cada padre en la generación $t$ con sus hijos en la generación $t+1$.
- Si un linaje no deja descendencia, su rama se corta (**extinción**).
- Conforme pasan las generaciones, la deriva genética reduce la diversidad hasta que **un único apellido conquista toda la población** (**fijación**).

Tomamos los **10 apellidos más frecuentes de Colombia** según el **Censo Nacional de Población y Vivienda (DANE 2018)**:
*Rodríguez, Martínez, García, Gómez, López, Hernández, Pérez, Sánchez, González y Ramírez*.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000031
md"**Modo de juego:** $(@bind modo_torneo Select([\"equitativo\" => \"🏁 Torneo Equitativo (1 representante por apellido, 10% prob. inicial c/u)\", \"dane\" => \"🇨🇴 Proporcional Censo DANE 2018 (Rodríguez arranca con ventaja poblacional)\"], default=\"equitativo\"))"

# ╔═╡ 00000001-0001-0001-0001-000000000049
md"**p\_herencia en el torneo** (0.5 = aleatorio, 1.0 = siempre hereda del origen): $(@bind p_her_torneo Slider(0.0:0.1:1.0, default=0.5, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000033
md"**Semilla de la partida (cámbiala para simular otra historia):** $(@bind semilla_torneo Slider(1:100, default=7, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000034
md"**Límite de generaciones a mostrar:** $(@bind max_gen_torneo Slider(10:5:40, default=25, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000035
begin
	# ── Datos oficiales Censo DANE 2018 (Top 10 Colombia) ──
	apellidos_colombia = [
		"Rodríguez", "Martínez", "García", "Gómez", "López",
		"Hernández", "Pérez", "Sánchez", "González", "Ramírez"
	]
	conteos_dane = [754991, 602593, 574856, 569534, 565394, 502626, 493764, 487112, 470505, 441179]
	pesos_dane = conteos_dane ./ sum(conteos_dane)
	
	colores_torneo = [
		"#1E88E5", # 1 Rodríguez (Azul)
		"#8E24AA", # 2 Martínez (Púrpura)
		"#FB8C00", # 3 García (Naranja)
		"#212121", # 4 Gómez (Negro)
		"#6D4C41", # 5 López (Marrón)
		"#E53935", # 6 Hernández (Rojo)
		"#43A047", # 7 Pérez (Verde)
		"#00ACC1", # 8 Sánchez (Turquesa)
		"#FDD835", # 9 González (Amarillo)
		"#3949AB"  # 10 Ramírez (Índigo)
	]
	
	N_torneo = 10
	rng_t = MersenneTwister(semilla_torneo)
	
	# Población inicial de la partida
	pob_init = if modo_torneo == "equitativo"
		collect(1:10)
	else
		sample(rng_t, 1:10, Weights(pesos_dane), N_torneo)
	end
	
	# ── Simular usando siguiente_generacion_grafo() del proyecto ──
	hist_genealogico = [copy(pob_init)]
	xs_conns = Float64[]
	ys_conns = Float64[]
	caidos = Dict{Int, Int}()
	
	for g in 1:max_gen_torneo
		pob_actual = hist_genealogico[g]
		
		# Usamos siguiente_generacion_grafo() con p_herencia del slider
		hijos = siguiente_generacion_grafo(pob_actual, rng_t; p_herencia=p_her_torneo)
		push!(hist_genealogico, hijos)
		
		# Para las líneas de parentesco, necesitamos una correspondencia visual.
		# siguiente_generacion_grafo produce N hijos, cada uno muestreado de la 
		# distribución completa. Dibujamos la conexión más probable (apellido del padre).
		for hijo_idx in 1:N_torneo
			ap_hijo = hijos[hijo_idx]
			# Buscar un padre con el mismo apellido para la conexión visual
			padres_posibles = findall(x -> x == ap_hijo, pob_actual)
			padre_idx = isempty(padres_posibles) ? rand(rng_t, 1:N_torneo) : rand(rng_t, padres_posibles)
			push!(xs_conns, padre_idx); push!(ys_conns, -(g - 1))
			push!(xs_conns, hijo_idx);  push!(ys_conns, -g)
			push!(xs_conns, NaN);       push!(ys_conns, NaN)
		end
		
		# Registrar extinciones
		sobrevivientes = unique(hijos)
		for ap in unique(hist_genealogico[1])
			if ap ∉ sobrevivientes && !haskey(caidos, ap)
				caidos[ap] = g
			end
		end
		
		if length(sobrevivientes) == 1
			break
		end
	end
	
	G_torneo = length(hist_genealogico) - 1
	ganador_id = hist_genealogico[end][1]
	hubo_fijacion = length(unique(hist_genealogico[end])) == 1
	
	# ── Dibujar la red genealógica ──
	p_arbol = plot(
		size = (800, max(460, 32 * (G_torneo + 1))),
		legend = :outerright,
		xlims = (0.3, N_torneo + 0.7),
		ylims = (-G_torneo - 0.7, 0.7),
		yticks = (-(0:G_torneo), ["Gen $g" for g in 0:G_torneo]),
		xticks = (1:N_torneo, ["Indiv. $i" for i in 1:N_torneo]),
		title = "Árbol Genealógico de Linajes (Wright-Fisher) — Colombia",
		grid = false,
		framestyle = :box
	)
	
	# Líneas de parentesco
	plot!(p_arbol, xs_conns, ys_conns, color = :black, alpha = 0.55, lw = 1.3, label = false)
	
	# Nodos agrupados por apellido
	xs_por_ap = [Float64[] for _ in 1:10]
	ys_por_ap = [Float64[] for _ in 1:10]
	
	for g in 0:G_torneo
		for i in 1:N_torneo
			ap = hist_genealogico[g + 1][i]
			push!(xs_por_ap[ap], i)
			push!(ys_por_ap[ap], -g)
		end
	end
	
	for ap in 1:10
		if !isempty(xs_por_ap[ap])
			scatter!(p_arbol, xs_por_ap[ap], ys_por_ap[ap],
				color = colores_torneo[ap],
				markersize = 10,
				markerstrokecolor = :black,
				markerstrokewidth = 1.1,
				label = "$(apellidos_colombia[ap])"
			)
		end
	end
	
	p_arbol
end

# ╔═╡ 00000001-0001-0001-0001-000000000036
begin
	ap_ganador_nom = apellidos_colombia[ganador_id]
	color_ganador = colores_torneo[ganador_id]
	
	extinciones_ordenadas = sort(collect(caidos), by = x -> x[2])
	filas_caidos = String[]
	for (pos, (ap_id, gen_ext)) in enumerate(extinciones_ordenadas)
		nom = apellidos_colombia[ap_id]
		push!(filas_caidos, "$pos. **$nom** (extinto en Gen $gen_ext)")
	end
	resumen_caidos = isempty(filas_caidos) ? "*Ninguno se extinguió aún.*" : join(filas_caidos, "<br>")
	
	md"""
	### 🏆 Resultados del Torneo
	
	| Categoría | Detalle |
	|:--|:--|
	| **Estado final** | $(hubo_fijacion ? "🎉 **¡FIJACIÓN ALCANZADA!**" : "⏳ No hubo fijación en $G_torneo generaciones") |
	| **Apellido Conquistador** | <span style="color:$(color_ganador); font-size:1.2em; font-weight:bold;">⬤ $(ap_ganador_nom)</span> |
	| **Generaciones hasta la victoria** | **$(G_torneo)** generaciones |
	| **p\_herencia utilizado** | **$(p_her_torneo)** |
	
	<br>
	
	**💀 Cronología de Extinción:**<br>
	$(resumen_caidos)
	
	> [!TIP]
	> **¿Por qué ocurre esto?** En una población finita sin mutación, el azar del muestreo (deriva genética) hace que los linajes compitan entre sí. Tarde o temprano, un único ancestro común acapara toda la población y todos los demás linajes se extinguen. ¡Prueba cambiando la semilla o el `p_herencia` para ver ganar a otro apellido!
	"""
end

# ╔═╡ 00000001-0001-0001-0001-000000000016
md"""
### 🎲 Múltiples Simulaciones: Distribución del Tiempo de Fijación

Corramos muchas simulaciones con `simular_fijacion()` para ver cuántas generaciones toma típicamente la fijación:
"""

# ╔═╡ 00000001-0001-0001-0001-000000000017
md"**Número de simulaciones**: $(@bind n_sims_slider Slider(50:50:1000, default=200, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000018
begin
	tiempos_fijacion = Int[]
	ganadores_hist = Int[]
	
	for s in 1:n_sims_slider
		rng_multi = MersenneTwister(s * 1000 + semilla_slider)
		
		# Construir población como en el script
		pob_multi = begin
			p = collect(1:num_apellidos)
			f = N_slider - num_apellidos
			if f > 0
				append!(p, rand(rng_multi, 1:num_apellidos, f))
			end
			shuffle(rng_multi, p)
		end
		
		# Llamar a simular_fijacion del proyecto
		res = simular_fijacion(pob_multi; parametros=parametros, rng=rng_multi)
		if res.fijacion
			push!(tiempos_fijacion, res.generaciones)
			push!(ganadores_hist, res.ganador)
		end
	end
	
	md"""
	**Resumen de $(n_sims_slider) simulaciones con `simular_fijacion()`:**
	- Fijaciones exitosas: **$(length(tiempos_fijacion))** / $(n_sims_slider)
	- Tiempo medio de fijación: **$(length(tiempos_fijacion) > 0 ? round(mean(tiempos_fijacion), digits=1) : "N/A")** generaciones
	- Tiempo mínimo: **$(length(tiempos_fijacion) > 0 ? minimum(tiempos_fijacion) : "N/A")**
	- Tiempo máximo: **$(length(tiempos_fijacion) > 0 ? maximum(tiempos_fijacion) : "N/A")**
	"""
end

# ╔═╡ 00000001-0001-0001-0001-000000000019
begin
	p3 = histogram(tiempos_fijacion,
		bins = max(10, div(length(tiempos_fijacion), 5)),
		xlabel = "Generaciones hasta fijación",
		ylabel = "Frecuencia",
		title = "Distribución del Tiempo de Fijación (N=$N_slider, $(num_apellidos) apellidos)",
		label = "Simulaciones", color = :steelblue, alpha = 0.7,
		size = (800, 400)
	)
	if length(tiempos_fijacion) > 0
		vline!(p3, [mean(tiempos_fijacion)], 
			label = "Media = $(round(mean(tiempos_fijacion), digits=1))",
			color = :red, linewidth = 2, linestyle = :dash)
	end
	p3
end

# ╔═╡ 00000001-0001-0001-0001-000000000020
begin
	if length(ganadores_hist) > 0
		conteo_ganadores = countmap(ganadores_hist)
		ap_ids = sort(collect(keys(conteo_ganadores)))
		victorias = [conteo_ganadores[a] for a in ap_ids]
		proporciones = victorias ./ sum(victorias)
		
		p4 = bar(
			["Ap. $a" for a in ap_ids], proporciones,
			xlabel = "Apellido", ylabel = "Proporción de victorias",
			title = "¿Qué apellido gana? ($(n_sims_slider) simulaciones)",
			label = false, color = palette(:tab10, length(ap_ids)),
			alpha = 0.8, size = (800, 400),
			ylims = (0, maximum(proporciones) * 1.2)
		)
		hline!(p4, [1/num_apellidos],
			label = "P uniforme (1/$(num_apellidos))",
			color = :red, linestyle = :dash, linewidth = 1.5)
		p4
	end
end

# ╔═╡ 00000001-0001-0001-0001-000000000021
md"""
### 🧮 Probabilidad Teórica de Fijación

En el modelo de Wright-Fisher neutral (sin selección), la probabilidad de que un apellido termine dominando es **proporcional a su frecuencia inicial**:

$$P(\text{fijación del apellido } i) = \frac{n_i}{N}$$

donde $n_i$ es el número de portadores iniciales y $N$ es el tamaño total.

El **tiempo esperado de fijación** escala como $\bar{T} \approx O(N)$ generaciones.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000022
begin
	conteo_inicial = countmap(estado_inicial)
	prob_teorica = Dict(ap => conteo_inicial[ap] / N_slider for ap in keys(conteo_inicial))
	
	lineas_prob = ["| Apellido | Portadores iniciales | Freq. inicial | P(fijación) teórica |",
	               "|:--|:--|:--|:--|"]
	for ap in sort(collect(keys(prob_teorica)))
		push!(lineas_prob, "| Ap. $ap | $(conteo_inicial[ap]) | $(round(prob_teorica[ap], digits=3)) | **$(round(prob_teorica[ap]*100, digits=1))%** |")
	end
	
	Markdown.parse(join(lineas_prob, "\n"))
end

# ╔═╡ 00000001-0001-0001-0001-000000000023
md"""
### 📝 Resumen

- El modelo usa las funciones de `src/stacional-model.jl` y `src/graphModel.jl`: `GrafoApellidos`, `siguiente_generacion_grafo()`, `ParametrosModelo` y `simular_fijacion()`
- La **topología del grafo de compatibilidad** define qué apellidos pueden interactuar: grafo completo (clásico), anillo (vecinos inmediatos), o mundo pequeño (anillo + atajos)
- El parámetro `p_herencia` controla la probabilidad de que el hijo herede el apellido del origen vs. el de la pareja
- La fijación es **inevitable** en poblaciones finitas: la diversidad sólo disminuye
- Poblaciones **más grandes** → más generaciones hasta fijación
- La **probabilidad de ganar** de cada apellido es proporcional a su **frecuencia inicial**
- Esto explica por qué ciertos apellidos comunes persisten, incluso sin ventaja reproductiva

---
*Notebook interactivo para el proyecto [Evolución y Extinción de Apellidos](https://github.com/spartan14799/Evolucion-y-Extincion-de-apellidos)*
"""

# ╔═╡ Cell order:
# ╟─00000001-0001-0001-0001-000000000001
# ╟─00000001-0001-0001-0001-000000000002
# ╟─00000001-0001-0001-0001-000000000003
# ╟─00000001-0001-0001-0001-000000000004
# ╟─00000001-0001-0001-0001-000000000005
# ╟─00000001-0001-0001-0001-000000000006
# ╟─00000001-0001-0001-0001-000000000007
# ╟─00000001-0001-0001-0001-000000000008
# ╟─00000001-0001-0001-0001-000000000040
# ╟─00000001-0001-0001-0001-000000000009
# ╟─00000001-0001-0001-0001-000000000010
# ╟─00000001-0001-0001-0001-000000000011
# ╟─00000001-0001-0001-0001-000000000012
# ╟─00000001-0001-0001-0001-000000000013
# ╟─00000001-0001-0001-0001-000000000014
# ╟─00000001-0001-0001-0001-000000000015
# ╟─00000001-0001-0001-0001-000000000041
# ╟─00000001-0001-0001-0001-000000000042
# ╟─00000001-0001-0001-0001-000000000043
# ╟─00000001-0001-0001-0001-000000000044
# ╟─00000001-0001-0001-0001-000000000045
# ╟─00000001-0001-0001-0001-000000000046
# ╟─00000001-0001-0001-0001-000000000047
# ╟─00000001-0001-0001-0001-000000000048
# ╟─00000001-0001-0001-0001-000000000030
# ╟─00000001-0001-0001-0001-000000000031
# ╟─00000001-0001-0001-0001-000000000049
# ╟─00000001-0001-0001-0001-000000000033
# ╟─00000001-0001-0001-0001-000000000034
# ╟─00000001-0001-0001-0001-000000000035
# ╟─00000001-0001-0001-0001-000000000036
# ╟─00000001-0001-0001-0001-000000000016
# ╟─00000001-0001-0001-0001-000000000017
# ╟─00000001-0001-0001-0001-000000000018
# ╟─00000001-0001-0001-0001-000000000019
# ╟─00000001-0001-0001-0001-000000000020
# ╟─00000001-0001-0001-0001-000000000021
# ╟─00000001-0001-0001-0001-000000000022
# ╟─00000001-0001-0001-0001-000000000023
