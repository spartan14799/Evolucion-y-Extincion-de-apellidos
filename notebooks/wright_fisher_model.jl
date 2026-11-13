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
	using Random, StatsBase, Plots, PlutoUI

	# ── Cargamos los módulos del proyecto tal cual están en src/ ──
	include(joinpath(@__DIR__, "..", "src", "stacional-model.jl"))
	# Esto nos da: ParametrosModelo, siguiente_generacion, simular_fijacion
end

# ╔═╡ 00000001-0001-0001-0001-000000000002
md"""
# 🏛️ Modelo Wright-Fisher para Apellidos
## Evolución y Fijación con Población Constante

---

### 📖 ¿Qué es este modelo?

El **modelo de Wright-Fisher** (originalmente diseñado para genética de poblaciones) modela cómo evoluciona la frecuencia de variantes en una **población de tamaño fijo**. En nuestro caso, los "alelos" son **apellidos**.

La idea central es (tal como la implementa [`src/stacional-model.jl`](../src/stacional-model.jl)):

> En cada generación, los `N` individuos se **mezclan al azar** (`shuffle`) y se agrupan en parejas de tamaño `k` (por defecto `k = 2`). De cada pareja sobrevive **un único apellido** (el del primer miembro tras la mezcla), que se propaga a todos los miembros del grupo. El tamaño de la población **nunca cambia**, pero la diversidad de apellidos **sólo puede disminuir**.

Esto es lo que hace la función `siguiente_generacion()` del proyecto.

El proceso avanza hacia la **fijación** de un solo apellido (estado absorbente), análogo a la **deriva genética**.

!!! info "Analogía biológica"
    Es el mismo fenómeno por el cual, sin selección natural, una mutación neutra termina fijándose o perdiéndose en una población finita. Aquí los "genes" son apellidos.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000003
md"""
### ⛓️ La Cadena de Markov

El estado del sistema en cada generación es el **vector de apellidos** de todos los individuos. Como la transición sólo depende del estado actual (mezcla aleatoria + selección de apellido ganador por pareja), el proceso es una **cadena de Markov**.

**Propiedades clave de `simular_fijacion()`:**
- **Espacio de estados**: Todas las posibles asignaciones de apellidos a `N` individuos
- **Transiciones**: `siguiente_generacion()` → shuffle → agrupar en parejas → propagar ganador
- **Estado absorbente**: Un único apellido en toda la población (`length(unique(poblacion)) == 1`)
- **Irreversibilidad**: La diversidad **nunca aumenta** → la fijación es **inevitable**
- **Invariante**: El tamaño de la población se conserva en cada paso (`N` constante)

La pregunta no es **si** un apellido domina, sino **cuándo** y **cuál**.
"""

# ╔═╡ 00000001-0001-0001-0001-000000000004
md"""
### 🎛️ Panel de Control

Ajusta los parámetros de la simulación (los mismos que usa `scripts/run_stacional_model.jl`):
"""

# ╔═╡ 00000001-0001-0001-0001-000000000005
md"**Población total** (N, debe ser par): $(@bind N_slider Slider(4:2:60, default=10, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000006
md"**Número de apellidos iniciales**: $(@bind num_apellidos_slider Slider(2:1:20, default=4, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000007
md"**Semilla aleatoria** (para reproducibilidad): $(@bind semilla_slider Slider(1:100, default=42, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000008
md"**Máximo de generaciones**: $(@bind max_gen_slider Slider(100:100:10000, default=5000, show_value=true))"

# ╔═╡ 00000001-0001-0001-0001-000000000009
begin
	# Ajustar apellidos al máximo posible (misma validación que el script)
	num_apellidos = min(num_apellidos_slider, N_slider)
	
	rng_sim = MersenneTwister(semilla_slider)

	# ── Construimos la población inicial ──
	# Esta lógica replica `configurar_poblacion()` de scripts/run_stacional_model.jl:
	# 1. Aseguramos que cada apellido (1:num_apellidos) aparezca al menos una vez
	# 2. El resto se rellena al azar
	# 3. Se mezcla
	estado_inicial = begin
		pob = collect(1:num_apellidos)
		faltantes = N_slider - num_apellidos
		if faltantes > 0
			append!(pob, rand(rng_sim, 1:num_apellidos, faltantes))
		end
		shuffle(rng_sim, pob)
	end
	
	# ── Usamos ParametrosModelo y simular_fijacion del proyecto ──
	parametros = ParametrosModelo(tamano_pareja = 2, max_generaciones = max_gen_slider)
	resultado = simular_fijacion(estado_inicial; parametros=parametros, rng=rng_sim)
	
	md"""
	### 📊 Resultado de la Simulación
	
	| Parámetro | Valor |
	|:--|:--|
	| Población (N) | **$(N_slider)** |
	| Apellidos iniciales | **$(num_apellidos)** |
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

El área total siempre es 1 (100% de la población se conserva — invariante de `siguiente_generacion()`).
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

Detalle de las primeras y últimas generaciones (mismo formato que `scripts/run_stacional_model.jl`):
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

# ╔═╡ 00000001-0001-0001-0001-000000000030
md"""
---

## 🎮 Minijuego: El Torneo de Linajes de Colombia
### Árbol Genealógico de Coalescencia (Censo DANE 2018)

Aquí simulamos el árbol genealógico completo generación por generación, exactamente como en el diagrama clásico de **Wright-Fisher**:
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

# ╔═╡ 00000001-0001-0001-0001-000000000032
md"**Dinámica de reproducción:** $(@bind dinamica_torneo Select([\"wf_clasico\" => \"🎲 Wright-Fisher Estándar (cada hijo elige padre al azar)\", \"parejas\" => \"👥 Modelo de Parejas del Proyecto (competencia en parejas k=2, stacional-model.jl)\"], default=\"wf_clasico\"))"

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
	
	hist_genealogico = [copy(pob_init)]
	xs_conns = Float64[]
	ys_conns = Float64[]
	caidos = Dict{Int, Int}()
	
	for g in 1:max_gen_torneo
		padres = Vector{Int}(undef, N_torneo)
		if dinamica_torneo == "wf_clasico"
			padres = rand(rng_t, 1:N_torneo, N_torneo)
		else
			orden = shuffle(rng_t, 1:N_torneo)
			for i in 1:2:N_torneo
				ganador_pareja = orden[i]
				padres[i] = ganador_pareja
				padres[i+1] = ganador_pareja
			end
		end
		
		hijos = [hist_genealogico[g][p] for p in padres]
		push!(hist_genealogico, hijos)
		
		# Trazar líneas padre -> hijo
		for (hijo_idx, padre_idx) in enumerate(padres)
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
	| **Dinámica utilizada** | $(dinamica_torneo == "wf_clasico" ? "Wright-Fisher estándar" : "Modelo de parejas (proyecto)") |
	
	<br>
	
	**💀 Cronología de Extinción:**<br>
	$(resumen_caidos)
	
	> [!TIP]
	> **¿Por qué ocurre esto?** En una población finita sin mutación, el azar del muestreo (deriva genética) hace que los linajes compitan entre sí. Tarde o temprano, un único ancestro común acapara toda la población y todos los demás linajes se extinguen. ¡Prueba cambiando la semilla para ver ganar a otro apellido!
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

- El modelo usa las funciones de `src/stacional-model.jl`: `ParametrosModelo`, `siguiente_generacion()` y `simular_fijacion()`
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
# ╟─00000001-0001-0001-0001-000000000009
# ╟─00000001-0001-0001-0001-000000000010
# ╟─00000001-0001-0001-0001-000000000011
# ╟─00000001-0001-0001-0001-000000000012
# ╟─00000001-0001-0001-0001-000000000013
# ╟─00000001-0001-0001-0001-000000000014
# ╟─00000001-0001-0001-0001-000000000015
# ╟─00000001-0001-0001-0001-000000000030
# ╟─00000001-0001-0001-0001-000000000031
# ╟─00000001-0001-0001-0001-000000000032
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
