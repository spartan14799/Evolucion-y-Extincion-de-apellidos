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

# ╔═╡ 00000002-0002-0002-0002-000000000001
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using Random, Distributions, Plots, PlutoUI, Printf

	# ── Cargamos los módulos del proyecto tal cual están en src/ ──
	include(joinpath(@__DIR__, "..", "src", "numericalMethods.jl"))
	# Esto nos da: biseccion, newton, secante

	include(joinpath(@__DIR__, "..", "src", "markovModels.jl"))
	# Esto nos da: ModeloPBGW, esperanza, es_extincion_segura, fgp,
	#              probabilidad_extincion, simular_generaciones, frecuencia_extincion
end

# ╔═╡ 00000002-0002-0002-0002-000000000002
md"""
# 🌳 Modelo de Bienaymé-Galton-Watson (PBGW)
## Extinción y Supervivencia de Apellidos como Proceso de Ramificación

---

### 📖 ¿Qué es este modelo?

El **proceso de Bienaymé-Galton-Watson** (BGW) es un modelo clásico de procesos de ramificación. Fue propuesto originalmente en el siglo XIX precisamente para estudiar nuestro problema: **la extinción de apellidos aristocráticos**.

La idea (implementada en [`src/markovModels.jl`](../src/markovModels.jl)):

> Un apellido comienza con `Z₀` portadores. Cada portador tiene un número **aleatorio** de hijos (que heredan el apellido), según una distribución de descendencia `{pₖ}`. La población de la siguiente generación es la **suma** de todos los hijos. Si llega a **cero**, el apellido se ha **extinguido** para siempre.

Formalmente, sea $Z_n$ el número de portadores en la generación $n$:

$$Z_{n+1} = \sum_{i=1}^{Z_n} X_i$$

donde cada $X_i$ es i.i.d. con distribución `pk` (el vector de `ModeloPBGW`).

Si $Z_n = 0$ → **extinción** (estado absorbente). Esto es lo que rastrea `simular_generaciones()`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000003
md"""
### ⛓️ La Cadena de Markov y la Función Generadora

El proceso BGW es una cadena de Markov sobre los enteros no negativos.

La herramienta fundamental es la **Función Generadora de Probabilidad** (FGP), implementada como `fgp(modelo, t)` en el proyecto:

$$g(t) = \mathbb{E}[t^X] = \sum_{k=0}^{\infty} p_k \, t^k, \quad t \in [0, 1]$$

**El resultado clave:** La **probabilidad de extinción** $q$ es la menor raíz de:

$$g(q) = q$$

Para encontrar $q$, la función `probabilidad_extincion()` define $h(t) = g(t) - t$ y busca la raíz usando `biseccion()` o `secante()` de `src/numericalMethods.jl`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000004
md"""
### 🔑 Los Tres Regímenes

El parámetro clave es la **media de descendencia** $m = \mathbb{E}[X]$, calculada por `esperanza(modelo)`:

| Régimen | Condición | Destino | $q$ | Función del proyecto |
|:--|:--|:--|:--|:--|
| 🔴 **Subcrítico** | $m < 1$ | Extinción **segura** | $q = 1$ | `es_extincion_segura() → true` |
| 🟡 **Crítico** | $m = 1$ | Extinción **segura** | $q = 1$ | `es_extincion_segura() → true` |
| 🟢 **Supercrítico** | $m > 1$ | Puede **sobrevivir** | $q < 1$ | `es_extincion_segura() → false` |

**Interpretación:** Si en promedio cada portador tiene menos de un hijo, el linaje se contrae y desaparece. Solo cuando $m > 1$ hay probabilidad positiva de pervivencia.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000005
md"""
### 🎛️ Panel de Control

Elige la distribución de natalidad y sus parámetros:
"""

# ╔═╡ 00000002-0002-0002-0002-000000000006
md"**Tipo de distribución de natalidad:** $(@bind dist_tipo Select([\"Poisson\" => \"Poisson (λ)\", \"Binomial\" => \"Binomial (n, p)\"], default=\"Poisson\"))"

# ╔═╡ 00000002-0002-0002-0002-000000000007
md"""
!!! info "¿Qué distribución elegir?"
    - **Poisson(λ):** Cada portador tiene en promedio λ hijos. No hay límite superior fijo. Es el modelo clásico (el que usa el constructor `ModeloPBGW(λ)` del proyecto). La media es $m = \lambda$.
    
    - **Binomial(n, p):** Cada portador tiene como máximo `n` hijos, cada uno con probabilidad `p`. La media es $m = n \cdot p$. Útil cuando hay un número máximo natural de hijos posibles. Se construye pasando el vector `pk` al constructor `ModeloPBGW(pk)`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000008
if dist_tipo == "Poisson"
	md"**Tasa de natalidad λ** (media de hijos): $(@bind lambda_slider Slider(0.1:0.1:4.0, default=1.5, show_value=true))"
else
	md""
end

# ╔═╡ 00000002-0002-0002-0002-000000000009
if dist_tipo == "Binomial"
	md"**Máximo de hijos (n):** $(@bind n_binom_slider Slider(1:1:15, default=5, show_value=true))"
else
	md""
end

# ╔═╡ 00000002-0002-0002-0002-000000000010
if dist_tipo == "Binomial"
	md"**Probabilidad por hijo (p):** $(@bind p_binom_slider Slider(0.05:0.05:1.0, default=0.4, show_value=true))"
else
	md""
end

# ╔═╡ 00000002-0002-0002-0002-000000000011
md"**Población inicial de portadores** ($Z_0$): $(@bind z0_slider Slider(1:1:20, default=1, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000012
md"**Generaciones a simular:** $(@bind gen_slider Slider(5:5:200, default=30, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000013
md"**Semilla aleatoria:** $(@bind seed_slider Slider(1:100, default=42, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000014
begin
	# ── Construir el modelo usando los constructores del proyecto ──
	if dist_tipo == "Poisson"
		# Usamos el constructor ModeloPBGW(tasa_natalidad) de markovModels.jl
		modelo = ModeloPBGW(lambda_slider; max_hijos=20)
		param_str = "Poisson(λ=$(lambda_slider))"
	else
		# Para Binomial, construimos el vector pk y usamos ModeloPBGW(pk)
		dist_binom = Distributions.Binomial(n_binom_slider, p_binom_slider)
		pk_binom = [pdf(dist_binom, k) for k in 0:n_binom_slider]
		pk_binom = pk_binom ./ sum(pk_binom)  # normalizar
		modelo = ModeloPBGW(pk_binom)
		param_str = "Binomial(n=$(n_binom_slider), p=$(p_binom_slider))"
	end
	
	# ── Usar las funciones del proyecto ──
	m_esperanza = esperanza(modelo)       # de markovModels.jl
	q_extincion = probabilidad_extincion(modelo)  # usa biseccion de numericalMethods.jl
	extincion_segura = es_extincion_segura(modelo) # de markovModels.jl
	
	regimen = if m_esperanza < 1.0
		"🔴 Subcrítico"
	elseif isapprox(m_esperanza, 1.0, atol=0.01)
		"🟡 Crítico"
	else
		"🟢 Supercrítico"
	end
	
	md"""
	### 📊 Análisis del Modelo: $(param_str)
	
	| Parámetro | Valor |
	|:--|:--|
	| Distribución | **$(param_str)** |
	| Media de descendencia (m) | **$(round(m_esperanza, digits=4))** — `esperanza(modelo)` |
	| Régimen | **$(regimen)** |
	| ¿Extinción segura? | **$(extincion_segura ? "Sí (m ≤ 1)" : "No (m > 1)")** — `es_extincion_segura(modelo)` |
	| Probabilidad de extinción (q) | **$(round(q_extincion, digits=6))** — `probabilidad_extincion(modelo)` |
	| Prob. de supervivencia (1-q) | **$(round(1 - q_extincion, digits=6))** |
	| Población inicial (Z₀) | **$(z0_slider)** |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000015
md"""
### 📉 Distribución de Descendencia $\{p_k\}$

El vector `modelo.pk` almacena las probabilidades. `pk[1]` = P(0 hijos), `pk[2]` = P(1 hijo), etc.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000016
begin
	# Mostrar la distribución pk del modelo
	pk_vec = modelo.pk
	ultimo_significativo = findlast(x -> x > 1e-6, pk_vec)
	k_rango = 0:(ultimo_significativo - 1)
	pk_plot = pk_vec[1:ultimo_significativo]
	
	p_dist = bar(
		collect(k_rango), pk_plot,
		xlabel = "Número de hijos (k)",
		ylabel = "Probabilidad P(X = k)",
		title = "Distribución de Descendencia — $(param_str)",
		label = false, color = :teal, alpha = 0.8,
		size = (800, 350), xticks = k_rango
	)
	vline!(p_dist, [m_esperanza], 
		label = "Media m = $(round(m_esperanza, digits=3))",
		color = :red, linewidth = 2, linestyle = :dash)
	p_dist
end

# ╔═╡ 00000002-0002-0002-0002-000000000017
md"""
### 📐 Función Generadora de Probabilidad g(t) — `fgp(modelo, t)`

La FGP es la herramienta con la que `probabilidad_extincion()` calcula $q$:
- Panel izquierdo: $g(t)$ vs la identidad $y = t$. El punto donde se cruzan (aparte de $t=1$) es $q$.
- Panel derecho: $h(t) = g(t) - t$. La raíz de $h$ es $q$, encontrada por `biseccion()`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000018
begin
	ts = range(0, 1, length=500)
	g_vals = [fgp(modelo, t) for t in ts]   # Usando fgp() del proyecto
	h_vals = g_vals .- ts
	
	p_fgp = plot(layout = (1, 2), size = (900, 400))
	
	# Panel izquierdo: g(t) vs identidad
	plot!(p_fgp[1], ts, g_vals, label = "g(t) = fgp(modelo, t)", color = :blue, linewidth = 2.5)
	plot!(p_fgp[1], ts, ts, label = "y = t", color = :gray, linewidth = 1.5, linestyle = :dash)
	if q_extincion < 1.0 - 1e-6
		scatter!(p_fgp[1], [q_extincion], [q_extincion],
			label = "q = $(round(q_extincion, digits=4))",
			color = :red, markersize = 8, markershape = :star5)
	end
	scatter!(p_fgp[1], [1.0], [1.0], label = "t = 1 (trivial)", color = :orange, markersize = 6)
	plot!(p_fgp[1], xlabel = "t", ylabel = "g(t)", title = "FGP y Punto Fijo", legend = :topleft)
	
	# Panel derecho: h(t) = g(t) - t
	plot!(p_fgp[2], ts, h_vals, label = "h(t) = g(t) - t", color = :purple, linewidth = 2.5)
	hline!(p_fgp[2], [0], label = false, color = :gray, linestyle = :dash)
	if q_extincion < 1.0 - 1e-6
		scatter!(p_fgp[2], [q_extincion], [0],
			label = "Raíz q = $(round(q_extincion, digits=4))",
			color = :red, markersize = 8, markershape = :star5)
	end
	plot!(p_fgp[2], xlabel = "t", ylabel = "h(t)", title = "Raíz de h(t) — biseccion()", legend = :topright)
	
	p_fgp
end

# ╔═╡ 00000002-0002-0002-0002-000000000019
md"""
### 🎲 Trayectorias del Linaje — `simular_generaciones()`

Cada trayectoria se genera llamando a `simular_generaciones(modelo, gen; poblacion_inicial=Z₀)` de `markovModels.jl`. Las rojas se extinguieron, las verdes sobrevivieron.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000020
md"**Número de trayectorias:** $(@bind n_traj_slider Slider(1:1:20, default=8, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000021
begin
	rng_traj = MersenneTwister(seed_slider)
	
	p_traj = plot(
		xlabel = "Generación", ylabel = "Portadores (Zₙ)",
		title = "Trayectorias con simular_generaciones() — $(param_str)",
		size = (850, 450), legend = :outerright,
		grid = true, gridstyle = :dash, gridalpha = 0.3
	)
	
	n_extintas = 0
	n_vivas = 0
	
	for i in 1:n_traj_slider
		# Llamamos directamente a simular_generaciones() del proyecto
		traj = simular_generaciones(modelo, gen_slider; 
			poblacion_inicial=z0_slider, tope_poblacion=10_000, rng=rng_traj)
		
		extinta = traj[end] == 0
		if extinta; n_extintas += 1; else; n_vivas += 1; end
		
		color = extinta ? :red : :green
		alpha_val = extinta ? 0.4 : 0.7
		
		plot!(p_traj, 0:gen_slider, traj,
			label = i <= 5 ? (extinta ? "Extinta #$i" : "Viva #$i") : false,
			color = color, alpha = alpha_val, linewidth = extinta ? 1.0 : 2.0)
	end
	
	p_traj
end

# ╔═╡ 00000002-0002-0002-0002-000000000022
md"""
**De $(n_traj_slider) trayectorias:** 🔴 $(n_extintas) extintas | 🟢 $(n_vivas) sobreviven

Probabilidad teórica: **q = $(round(q_extincion, digits=4))** → frecuencia observada: **$(round(n_extintas / n_traj_slider, digits=4))**
"""

# ╔═╡ 00000002-0002-0002-0002-000000000023
md"""
### 🎯 Validación Monte Carlo — `frecuencia_extincion()`

Comparemos la `probabilidad_extincion()` (teórica, bisección) contra `frecuencia_extincion()` (empírica, Monte Carlo). Ambas funciones están en `markovModels.jl`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000024
md"**Simulaciones Monte Carlo:** $(@bind n_mc_slider Slider(500:500:10000, default=3000, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000025
begin
	rng_mc = MersenneTwister(seed_slider + 7777)
	
	# Llamamos a frecuencia_extincion() del proyecto
	freq_mc = frecuencia_extincion(modelo, gen_slider;
		n_sim=n_mc_slider, poblacion_inicial=z0_slider,
		tope_poblacion=10_000, rng=rng_mc)
	
	diferencia = abs(freq_mc - q_extincion)
	
	md"""
	| Medida | Valor | Función |
	|:--|:--|:--|
	| q teórica | **$(round(q_extincion, digits=6))** | `probabilidad_extincion(modelo)` |
	| q empírica (n=$(n_mc_slider)) | **$(round(freq_mc, digits=6))** | `frecuencia_extincion(modelo, $(gen_slider))` |
	| Diferencia | **$(round(diferencia, digits=6))** | |
	| Concordancia | **$(diferencia < 0.05 ? "✅ Excelente" : (diferencia < 0.1 ? "⚠️ Aceptable" : "❌ Divergente"))** | |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000026
md"""
### 🔬 Transición de Fase: Extinción vs. Supervivencia

¿Cómo cambia $q$ al variar la media $m$? En $m = 1$ ocurre una **transición de fase** abrupta.

Cada punto se calcula con `probabilidad_extincion(ModeloPBGW(λ))`:
"""

# ╔═╡ 00000002-0002-0002-0002-000000000027
begin
	lambdas = range(0.1, 4.0, length=100)
	qs = Float64[]
	for lam in lambdas
		# Construimos modelos Poisson y calculamos q con el proyecto
		mod_temp = ModeloPBGW(lam; max_hijos=20)
		push!(qs, probabilidad_extincion(mod_temp))
	end
	
	p_phase = plot(
		lambdas, qs,
		xlabel = "Tasa de natalidad (λ = media m)",
		ylabel = "Probabilidad de extinción (q)",
		title = "Transición de Fase: q vs m",
		label = "q(m) — probabilidad_extincion()",
		color = :darkblue, linewidth = 3,
		size = (800, 400), ylims = (-0.05, 1.1), grid = true
	)
	vline!(p_phase, [1.0], label = "m = 1 (umbral)", color = :red, linestyle = :dash, linewidth = 2)
	
	if dist_tipo == "Poisson"
		scatter!(p_phase, [lambda_slider], [q_extincion],
			label = "Modelo actual (λ=$(lambda_slider))",
			color = :orange, markersize = 10, markershape = :star5)
	end
	
	annotate!(p_phase, 0.5, 0.3, text("Subcrítico\n(extinción segura)", 9, :red))
	annotate!(p_phase, 2.5, 0.3, text("Supercrítico\n(supervivencia posible)", 9, :green))
	
	p_phase
end

# ╔═╡ 00000002-0002-0002-0002-000000000028
md"""
### 📊 Tabla Resumen para Distintas Tasas de Natalidad
"""

# ╔═╡ 00000002-0002-0002-0002-000000000029
begin
	lineas_tabla = [
		"| λ (media m) | Régimen | q (extinción) | 1-q (supervivencia) |",
		"|:--|:--|:--|:--|"
	]
	for lam in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
		mod_t = ModeloPBGW(lam; max_hijos=20)
		m_t = esperanza(mod_t)
		q_t = probabilidad_extincion(mod_t)
		reg_t = m_t < 1.0 ? "🔴 Subcrítico" : (isapprox(m_t, 1.0, atol=0.01) ? "🟡 Crítico" : "🟢 Supercrítico")
		push!(lineas_tabla, "| $(lam) | $(reg_t) | $(round(q_t, digits=4)) | $(round(1-q_t, digits=4)) |")
	end
	Markdown.parse(join(lineas_tabla, "\n"))
end

# ╔═╡ 00000002-0002-0002-0002-000000000030
md"""
### 🧮 Probabilidad de Absorción con Z₀ > 1

Cuando el linaje empieza con $Z_0 > 1$ portadores, cada uno genera un sublinaje **independiente**. Todos deben extinguirse para que el apellido desaparezca:

$$q_{Z_0} = q^{Z_0}$$

Esto se puede verificar con `frecuencia_extincion(modelo, gen; poblacion_inicial=Z₀)`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000031
begin
	z0_vals = 1:10
	q_z0_vals = [q_extincion^z for z in z0_vals]
	
	p_absorcion = bar(
		collect(z0_vals), q_z0_vals,
		xlabel = "Población inicial Z₀",
		ylabel = "P(extinción) = q^(Z₀)",
		title = "Efecto de Z₀ — $(param_str) (q=$(round(q_extincion, digits=3)))",
		label = false,
		color = [q < 0.5 ? :green : :red for q in q_z0_vals],
		alpha = 0.8, size = (800, 400), xticks = z0_vals, ylims = (0, 1.05)
	)
	hline!(p_absorcion, [0.5], label = "50%", color = :gray, linestyle = :dash)
	for (z, q_z) in zip(z0_vals, q_z0_vals)
		annotate!(p_absorcion, z, q_z + 0.03, text("$(round(q_z*100, digits=1))%", 7, :center))
	end
	p_absorcion
end

# ╔═╡ 00000002-0002-0002-0002-000000000032
md"""
### 📝 Resumen y Conclusiones

**Funciones del proyecto utilizadas en este notebook:**

| Función | Archivo | Propósito |
|:--|:--|:--|
| `ModeloPBGW(λ)` / `ModeloPBGW(pk)` | `src/markovModels.jl` | Construir modelo con dist. Poisson o vector explícito |
| `esperanza(modelo)` | `src/markovModels.jl` | Media de descendencia $m$ |
| `es_extincion_segura(modelo)` | `src/markovModels.jl` | Verificar si $m \leq 1$ |
| `fgp(modelo, t)` | `src/markovModels.jl` | Función generadora $g(t)$ |
| `probabilidad_extincion(modelo)` | `src/markovModels.jl` | Resolver $g(t) = t$ con bisección |
| `simular_generaciones(modelo, gen)` | `src/markovModels.jl` | Simular una trayectoria |
| `frecuencia_extincion(modelo, gen)` | `src/markovModels.jl` | Monte Carlo: fracción de extintos |
| `biseccion(h, a, b)` | `src/numericalMethods.jl` | Búsqueda de raíces (usado internamente) |
| `secante(h, a, b)` | `src/numericalMethods.jl` | Búsqueda de raíces (alternativa) |

**Conclusiones clave:**
- Todo depende de $m$: si $m \leq 1$, extinción segura; si $m > 1$, supervivencia posible
- La FGP + bisección nos dan la probabilidad exacta de extinción
- Monte Carlo confirma la teoría: `frecuencia_extincion ≈ probabilidad_extincion`
- Con $Z_0 > 1$ portadores iniciales, la extinción baja exponencialmente: $q^{Z_0}$

---
*Notebook interactivo para el proyecto [Evolución y Extinción de Apellidos](https://github.com/spartan14799/Evolucion-y-Extincion-de-apellidos)*
"""

# ╔═╡ Cell order:
# ╟─00000002-0002-0002-0002-000000000001
# ╟─00000002-0002-0002-0002-000000000002
# ╟─00000002-0002-0002-0002-000000000003
# ╟─00000002-0002-0002-0002-000000000004
# ╟─00000002-0002-0002-0002-000000000005
# ╟─00000002-0002-0002-0002-000000000006
# ╟─00000002-0002-0002-0002-000000000007
# ╟─00000002-0002-0002-0002-000000000008
# ╟─00000002-0002-0002-0002-000000000009
# ╟─00000002-0002-0002-0002-000000000010
# ╟─00000002-0002-0002-0002-000000000011
# ╟─00000002-0002-0002-0002-000000000012
# ╟─00000002-0002-0002-0002-000000000013
# ╠═00000002-0002-0002-0002-000000000014
# ╟─00000002-0002-0002-0002-000000000015
# ╠═00000002-0002-0002-0002-000000000016
# ╟─00000002-0002-0002-0002-000000000017
# ╠═00000002-0002-0002-0002-000000000018
# ╟─00000002-0002-0002-0002-000000000019
# ╟─00000002-0002-0002-0002-000000000020
# ╠═00000002-0002-0002-0002-000000000021
# ╟─00000002-0002-0002-0002-000000000022
# ╟─00000002-0002-0002-0002-000000000023
# ╟─00000002-0002-0002-0002-000000000024
# ╠═00000002-0002-0002-0002-000000000025
# ╟─00000002-0002-0002-0002-000000000026
# ╠═00000002-0002-0002-0002-000000000027
# ╟─00000002-0002-0002-0002-000000000028
# ╠═00000002-0002-0002-0002-000000000029
# ╟─00000002-0002-0002-0002-000000000030
# ╠═00000002-0002-0002-0002-000000000031
# ╟─00000002-0002-0002-0002-000000000032
