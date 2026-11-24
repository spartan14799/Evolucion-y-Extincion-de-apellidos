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
	include(joinpath(@__DIR__, "..", "src", "core", "numericalMethods.jl"))
	# Esto nos da: biseccion, newton, secante

	include(joinpath(@__DIR__, "..", "src", "bgwModel.jl"))
	# Esto nos da: ModeloPBGW, esperanza, es_extincion_segura, fgp,
	#              probabilidad_extincion, simular_generaciones, frecuencia_extincion
end

# ╔═╡ 00000002-0002-0002-0002-000000000002
md"""
# Modelo de Bienaymé-Galton-Watson (PBGW)
## Extinción y Supervivencia de Apellidos como Proceso de Ramificación

---

### ¿Qué es este modelo?

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
### La Cadena de Markov y la Función Generadora

El proceso BGW es una cadena de Markov sobre los enteros no negativos.

La herramienta fundamental es la **Función Generadora de Probabilidad** (FGP), implementada como `fgp(modelo, t)` en el proyecto:

$$g(t) = \mathbb{E}[t^X] = \sum_{k=0}^{\infty} p_k \, t^k, \quad t \in [0, 1]$$

**El resultado clave:** La **probabilidad de extinción** $q$ es la menor raíz de:

$$g(q) = q$$

Para encontrar $q$, la función `probabilidad_extincion()` define $h(t) = g(t) - t$ y busca la raíz usando `biseccion()` o `secante()` de `src/numericalMethods.jl`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000004
md"""
### Los Tres Regímenes

El parámetro clave es la **media de descendencia** $m = \mathbb{E}[X]$, calculada por `esperanza(modelo)`:

| Régimen | Condición | Destino | $q$ | Función del proyecto |
|:--|:--|:--|:--|:--|
| **Subcrítico** | $m < 1$ | Extinción **segura** | $q = 1$ | `es_extincion_segura() → true` |
| **Crítico** | $m = 1$ | Extinción **segura** | $q = 1$ | `es_extincion_segura() → true` |
| **Supercrítico** | $m > 1$ | Puede **sobrevivir** | $q < 1$ | `es_extincion_segura() → false` |

**Interpretación:** Si en promedio cada portador tiene menos de un hijo, el linaje se contrae y desaparece. Solo cuando $m > 1$ hay probabilidad positiva de pervivencia.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000005
md"""
### Panel de Control

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
		"Subcrítico"
	elseif isapprox(m_esperanza, 1.0, atol=0.01)
		"Crítico"
	else
		"Supercrítico"
	end
	
	md"""
	### Análisis del Modelo: $(param_str)
	
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
### Distribución de Descendencia $\{p_k\}$

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
### Función Generadora de Probabilidad g(t) — `fgp(modelo, t)`

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
### Trayectorias del Linaje — `simular_generaciones()`

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
**De $(n_traj_slider) trayectorias:** $(n_extintas) extintas | $(n_vivas) sobreviven

Probabilidad teórica: **q = $(round(q_extincion, digits=4))** → frecuencia observada: **$(round(n_extintas / n_traj_slider, digits=4))**
"""

# ╔═╡ 00000002-0002-0002-0002-000000000023
md"""
### Validación Monte Carlo — `frecuencia_extincion()`

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
	| Concordancia | **$(diferencia < 0.05 ? "Excelente" : (diferencia < 0.1 ? "Aceptable" : "Divergente"))** | |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000026
md"""
### Transición de Fase: Extinción vs. Supervivencia

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
### Tabla Resumen para Distintas Tasas de Natalidad
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
		reg_t = m_t < 1.0 ? "Subcrítico" : (isapprox(m_t, 1.0, atol=0.01) ? "Crítico" : "Supercrítico")
		push!(lineas_tabla, "| $(lam) | $(reg_t) | $(round(q_t, digits=4)) | $(round(1-q_t, digits=4)) |")
	end
	Markdown.parse(join(lineas_tabla, "\n"))
end

# ╔═╡ 00000002-0002-0002-0002-000000000030
md"""
### Probabilidad de Absorción con Z₀ > 1

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
### Resumen y Conclusiones

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

# ╔═╡ 00000002-0002-0002-0002-000000000040
md"""
---

# Simulador: ¿Cuándo se extinguirán los apellidos de Colombia?
## Proceso de Bienaymé-Galton-Watson con Datos Oficiales del Censo DANE 2018

### Contexto demográfico real:
En Colombia, según proyecciones y datos demográficos recientes (DANE y CEPAL), la tasa global de fecundidad ha caído a cerca de **1.2 - 1.5 hijos por mujer** (en varias regiones urbanas ya se sitúa en torno a 1.0).
Dado que tradicionalmente el primer apellido lo transmiten los hijos varones:
- Si una mujer tiene en promedio 1 hijo, aproximadamente **0.5 son varones**.
- Por lo tanto, cada hombre que porta el apellido engendra en promedio **$\lambda \approx 0.5$ hijos varones** que continuarán su linaje.
- Dado que $m = \lambda = 0.5 \le 1$, el linaje se encuentra en **régimen subcrítico**: **la extinción matemática es 100% segura ($q = 1$)**.

La pregunta fundamental ya no es *si* se va a extinguir, sino: **¿en cuántas generaciones ocurrirá la extinción para los apellidos más populosos de Colombia?**
"""

# ╔═╡ 00000002-0002-0002-0002-000000000041
md"**Selecciona un apellido de Colombia:** $(@bind ap_col_select Select([\"Rodríguez\" => \"1. Rodríguez (754.991 personas)\", \"Martínez\" => \"2. Martínez (602.593 personas)\", \"García\" => \"3. García (574.856 personas)\", \"Gómez\" => \"4. Gómez (569.534 personas)\", \"López\" => \"5. López (565.394 personas)\", \"Hernández\" => \"6. Hernández (502.626 personas)\", \"Pérez\" => \"7. Pérez (493.764 personas)\", \"Sánchez\" => \"8. Sánchez (487.112 personas)\", \"González\" => \"9. González (470.505 personas)\", \"Ramírez\" => \"10. Ramírez (441.179 personas)\"], default=\"Rodríguez\"))"

# ╔═╡ 00000002-0002-0002-0002-000000000042
md"**Tasa de natalidad masculina $\lambda$** (hijos varones promedio por hombre): $(@bind lambda_col_slider Slider(0.1:0.05:1.2, default=0.5, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000043
md"**Evaluar probabilidad de extinción acumulada a las generaciones ($X$):** $(@bind gen_eval_slider Slider(5:1:40, default=20, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000044
md"**Simulaciones Monte Carlo:** $(@bind n_sims_col_slider Slider(50:25:500, default=100, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000045
md"**Semilla aleatoria:** $(@bind semilla_col_slider Slider(1:100, default=42, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000046
begin
	apellidos_colombia = [
		"Rodríguez", "Martínez", "García", "Gómez", "López",
		"Hernández", "Pérez", "Sánchez", "González", "Ramírez"
	]
	conteos_dane = [754991, 602593, 574856, 569534, 565394, 502626, 493764, 487112, 470505, 441179]
	dane_dict = Dict(zip(apellidos_colombia, conteos_dane))
	
	z0_seleccionado = dane_dict[ap_col_select]
	
	# Instanciamos el modelo de Bienaymé-Galton-Watson llamando a ModeloPBGW de src/markovModels.jl
	modelo_colombia = ModeloPBGW(lambda_col_slider; max_hijos = 15)
	
	# Tabla comparativa resaltando el apellido activo (sin emojis)
	filas_tabla = [
		"| Puesto | Apellido | Portadores DANE 2018 (\$Z_0\$) | Estado en Simulación |",
		"|:---|:---|:---|:---|"
	]
	for (idx, (ap, cnt)) in enumerate(zip(apellidos_colombia, conteos_dane))
		marca = ap == ap_col_select ? "**SELECCIONADO**" : "En espera"
		ap_fmt = ap == ap_col_select ? "**$ap**" : ap
		push!(filas_tabla, "| $idx | $ap_fmt | $(cnt) | $marca |")
	end
	Markdown.parse(join(filas_tabla, "\n"))
end

# ╔═╡ 00000002-0002-0002-0002-000000000047
md"""
### Trucos Numéricos y Tipo de Monte Carlo Utilizado

#### 1. Implementación Canónica del Código Fuente vs. Propiedad Reproductiva de Poisson
En el código fuente (`src/markovModels.jl`), la función `simular_generaciones` itera individuo por individuo sorteando su descendencia a través de la distribución categórica discreta `Categorical(modelo.pk)`:
```julia
# Enfoque canónico del código fuente (src/markovModels.jl):
for _ in 1:poblacion
    hijos_totales += rand(rng, dist_hijos) - 1
end
```
**El truco matemático para optimización analítica:**
Cuando la descendencia individual sigue una distribución de Poisson con media $\lambda$ ($X_i \stackrel{iid}{\sim} \text{Poisson}(\lambda)$), aprovechamos la **propiedad reproductiva de Poisson**: la suma de variables Poisson independientes es idéntica en distribución a una Poisson con la suma de las medias:
$$\sum_{i=1}^{Z_n} X_i \sim \text{Poisson}(Z_n \cdot \lambda)$$
Esto permite que para poblaciones grandes como las del Censo DANE ($Z_0 \approx 750.000$), la distribución exacta pueda generarse tanto sumando los individuos como evaluando directamente en tiempo $O(1)$. Ambas formulaciones son matemáticamente equivalentes y convergen a los mismos estimadores estadísticos.

#### 2. ¿Qué tipo de Monte Carlo usamos? (Forward MC vs MCMC)
- **Monte Carlo Directo / Hacia Adelante (Forward Monte Carlo):** Es el método empleado en este módulo. Simulamos realizaciones estocásticas directas del proceso temporal $\{Z_n\}_{n \ge 0}$ a lo largo de las generaciones hacia el futuro partiendo del estado conocido $Z_0$.
- **¿Por qué NO usamos MCMC (Markov Chain Monte Carlo - Metropolis/Gibbs)?**
  MCMC se utiliza para muestrear de una distribución estacionaria $\pi$ compleja o una distribución a posteriori bayesiana. En un proceso de ramificación subcrítico ($\lambda < 1$), los estados $\{1, 2, \dots\}$ son **transitorios** y el único estado recurrente es el $0$ (absorbente). No hay una distribución estacionaria no nula; el interés radica en la **trayectoria dinámica de extinción**, la cual se modela mediante Forward Monte Carlo con `simular_generaciones` de `src/markovModels.jl`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000048
begin
	# Simulación Monte Carlo llamando directamente a las funciones de src/markovModels.jl
	function correr_montecarlo_colombia(modelo::ModeloPBGW, z0::Int, n_sims::Int, gen_eval::Int, semilla::Int)
		rng = MersenneTwister(semilla)
		max_gen_sim = max(gen_eval + 15, 35)
		tiempos_ext = Int[]
		trayectorias_muestra = Vector{Vector{Int}}()
		extintos_en_X = 0
		
		for s in 1:n_sims
			# Llamada directa al método fuente simular_generaciones
			traj = simular_generaciones(modelo, max_gen_sim; poblacion_inicial = z0, rng = rng)
			if s <= 6
				push!(trayectorias_muestra, copy(traj))
			end
			idx_cero = findfirst(==(0), traj)
			if idx_cero !== nothing
				gen_ext = idx_cero - 1
				push!(tiempos_ext, gen_ext)
				if gen_ext <= gen_eval
					extintos_en_X += 1
				end
			end
		end
		
		prob_ext_X = extintos_en_X / n_sims
		mediana_ext = isempty(tiempos_ext) ? 0 : round(Int, median(tiempos_ext))
		media_ext = isempty(tiempos_ext) ? 0.0 : mean(tiempos_ext)
		
		return (
			tiempos_ext = tiempos_ext,
			trayectorias_muestra = trayectorias_muestra,
			extintos_en_X = extintos_en_X,
			prob_ext_X = prob_ext_X,
			mediana_ext = mediana_ext,
			media_ext = media_ext,
			max_gen_sim = max_gen_sim
		)
	end

	res_sim_col = correr_montecarlo_colombia(modelo_colombia, z0_seleccionado, n_sims_col_slider, gen_eval_slider, semilla_col_slider)
	tiempos_ext = res_sim_col.tiempos_ext
	trayectorias_muestra = res_sim_col.trayectorias_muestra
	extintos_en_X = res_sim_col.extintos_en_X
	prob_ext_X = res_sim_col.prob_ext_X
	mediana_ext = res_sim_col.mediana_ext
	media_ext = res_sim_col.media_ext
	max_gen_sim = res_sim_col.max_gen_sim
	
	# Validación teórica llamando directamente a src/markovModels.jl
	m_esp = esperanza(modelo_colombia)
	es_segura = es_extincion_segura(modelo_colombia)
	q_teorico = probabilidad_extincion(modelo_colombia)
	
	# Estimación teórica continua aproximada: T ≈ ln(Z0) / ln(1/m)
	t_teorico = m_esp < 1.0 ? log(z0_seleccionado) / log(1.0 / m_esp) : Inf
	
	md"""
	### Resultados del Monte Carlo para **$(ap_col_select)**
	
	| Métrica | Valor Estimado | Significado y Validación Teórica (markovModels.jl) |
	|:---|:---|:---|
	| **Población Inicial ($Z_0$)** | **$(z0_seleccionado)** | Portadores según Censo DANE 2018 |
	| **Esperanza de descendencia ($m = E[X]$)** | **$(round(m_esp, digits=3))** | Evaluado con `esperanza(modelo)` |
	| **Régimen del Proceso** | **$(es_segura ? "Subcrítico / Crítico (m ≤ 1.0)" : "Supercrítico (m > 1.0)")** | Evaluado con `es_extincion_segura(modelo)` |
	| **Probabilidad teórica asintótica ($q$)** | **$(round(q_teorico * 100, digits=2))%** | Calculado con `probabilidad_extincion(modelo)` |
	| **Probabilidad de extinción en $\le $(gen_eval_slider) generaciones** | <span style="font-size:1.3em; font-weight:bold; color:$(prob_ext_X > 0.5 ? "crimson" : "steelblue");">$(round(prob_ext_X * 100, digits=2))%</span> | Fracción observada en Monte Carlo con `simular_generaciones` |
	| **Generación Mediana de Extinción** | **Gen $(mediana_ext)** (~$(mediana_ext * 28) años) | El 50% de las realizaciones estocásticas se extinguieron antes de este año |
	| **Generación Media de Extinción** | **Gen $(round(media_ext, digits=1))** | Promedio aritmético de las simulaciones |
	| **Aproximación Teórica Determinista** | **Gen $(round(t_teorico, digits=1))** | Fórmula de decaimiento geométrico: $T \approx \frac{\ln(Z_0)}{\ln(1/m)}$ |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000049
begin
	# Histograma de Tiempos de Extinción
	if !isempty(tiempos_ext)
		p_hist_ext = histogram(
			tiempos_ext,
			bins = range(minimum(tiempos_ext)-0.5, maximum(tiempos_ext)+0.5, step=1),
			xlabel = "Generación de Extinción (1 gen ≈ 25 - 30 años)",
			ylabel = "Frecuencia de Simulaciones",
			title = "Distribución del Tiempo de Extinción para $(ap_col_select) (n=$(n_sims_col_slider))",
			color = "#e53935", alpha = 0.75, legend = :topright,
			size = (800, 380)
		)
		vline!(p_hist_ext, [mediana_ext], label = "Mediana = Gen $(mediana_ext)", color = :black, linewidth = 2.5, linestyle = :dash)
		vline!(p_hist_ext, [gen_eval_slider], label = "Umbral evaluado X = Gen $(gen_eval_slider) ($(round(prob_ext_X*100, digits=1))% extinto)", color = :blue, linewidth = 2)
		p_hist_ext
	else
		plot(title="No hubo extinciones (λ >= 1.0)")
	end
end

# ╔═╡ 00000002-0002-0002-0002-000000000050
begin
	# Trayectorias en Escala Logarítmica
	p_log = plot(
		xlabel = "Generación", ylabel = "Portadores (escala log10)",
		title = "Decaimiento Estocástico de Linajes — $(ap_col_select) ($(@sprintf("%.2f", lambda_col_slider)) varones/varón)",
		yaxis = :log10, size = (800, 400), legend = :bottomleft,
		grid = true, gridstyle = :dash, gridalpha = 0.3
	)
	
	# Curva esperada E[Z_n] = Z0 * lambda^n
	gens_curva = 0:min(max_gen_sim, 35)
	curva_esperada = [max(1.0, z0_seleccionado * (lambda_col_slider^g)) for g in gens_curva]
	plot!(p_log, gens_curva, curva_esperada, label = "Media esperada E[Z_n] = Z₀ · λⁿ", color = :black, linewidth = 2.5, linestyle = :dash)
	
	for (idx, tr) in enumerate(trayectorias_muestra)
		tr_plot = [max(0.5, Float64(x)) for x in tr]
		plot!(p_log, 0:(length(tr_plot)-1), tr_plot, label = idx == 1 ? "Simulaciones individuales" : false, alpha = 0.7, linewidth = 1.5)
	end
	hline!(p_log, [1.0], label = "1 portador (umbral crítico)", color = :red, linestyle = :dot, linewidth = 1.5)
	p_log
end

# ╔═╡ 00000002-0002-0002-0002-000000000051
begin
	# Comparativa de los 10 apellidos llamando directamente a simular_generaciones de src/markovModels.jl
	function calcular_medianas_10(modelo::ModeloPBGW)
		tiempos_medios = Float64[]
		rng_comp = MersenneTwister(123)
		for cnt in conteos_dane
			tiempos_i = Int[]
			for _ in 1:30
				traj_i = simular_generaciones(modelo, 35; poblacion_inicial = cnt, rng = rng_comp)
				idx_0 = findfirst(==(0), traj_i)
				push!(tiempos_i, idx_0 !== nothing ? idx_0 - 1 : 36)
			end
			push!(tiempos_medios, median(tiempos_i))
		end
		return tiempos_medios
	end

	tiempos_medios_10 = calcular_medianas_10(modelo_colombia)
	colores_barras = [ap == ap_col_select ? "#e53935" : "#1e88e5" for ap in apellidos_colombia]
	
	# Orden inverso para que Rodríguez (Puesto 1) aparezca arriba y se lean bien los nombres
	indices_bar = 1:10
	aps_rev = reverse(apellidos_colombia)
	tiempos_rev = reverse(tiempos_medios_10)
	colores_rev = reverse(colores_barras)

	p_comp = bar(
		indices_bar, tiempos_rev,
		orientation = :horizontal,
		yticks = (indices_bar, aps_rev),
		xlabel = "Generaciones estimadas hasta extinción",
		ylabel = "Apellido",
		title = "¿Cuánto tarda en extinguirse cada apellido? (λ = $(@sprintf("%.2f", lambda_col_slider)))",
		color = colores_rev, alpha = 0.85, legend = false,
		size = (800, 420), xlims = (0, maximum(tiempos_medios_10) * 1.35)
	)
	
	for (idx, t_med) in enumerate(tiempos_rev)
		annotate!(p_comp, t_med + 0.4, idx, text("Gen $(round(Int, t_med)) (~$(round(Int, t_med*28)) años)", 8, :left, :black, Sys.iswindows() ? "arialbd" : "sans-serif"))
	end
	p_comp
end

# ╔═╡ 00000002-0002-0002-0002-000000000052
md"""
> [!TIP]
> **La paradoja de los apellidos masivos:**
> Observa que **Rodríguez** (754.991 personas) tarda aproximadamente **20 generaciones** en extinguirse, mientras que **Ramírez** (441.179 personas, con más de 300.000 portadores menos) se extingue en **19 generaciones**.
> 
> **¿Por qué solo hay 1 generación de diferencia con 300.000 personas de ventaja?**
> Porque en un proceso de ramificación subcrítico, la población decae **geométricamente** ($Z_n \approx Z_0 \lambda^n$). Al despejar el tiempo de extinción:
> $$T_{\text{ext}} \approx \frac{\ln(Z_0)}{\ln(1/\lambda)}$$
> El tiempo crece con el **logaritmo** de la población inicial: ¡duplicar la cantidad de portadores apenas le otorga al apellido **una generación más de vida**!
"""

# ╔═╡ 00000002-0002-0002-0002-000000000060
md"""
---

# La Cadena de Markov por Dentro: Caracterización Explícita
## Estados, Matriz de Transición P y Flujo de Probabilidad

En las secciones anteriores hemos analizado las estadísticas globales del modelo. Aquí abrimos la máquina para ver el motor formal: **la cadena de Markov en tiempo discreto** subyacente.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000061
md"""
### 1. Definición del Espacio de Estados
En el modelo matemático riguroso, el espacio de estados es infinito numerable: $S = \{0, 1, 2, 3, \ldots\}$, donde cada estado representa la cantidad exacta de portadores vivos del apellido.

Para visualizar y calcular la matriz de forma explícita e interactiva, consideramos el **espacio de estados truncado** $S_M = \{0, 1, 2, \ldots, M\}$:
- **Estado 0 (Extinción):** Es un **estado absorbente**. Si el linaje llega a 0 personas, nadie puede reproducirse; la probabilidad de permanecer en 0 es 1 ($P_{0,0} = 1$).
- **Estados $1, 2, \dots, M-1$:** Son **estados transitorios**. El linaje tiene individuos activos, pero eventualmente la masa de probabilidad los abandonará.
- **Estado $M$ ($M+$ individuos):** Estado de acumulación superior que agrupa a todos los linajes con $M$ o más portadores.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000062
md"**Tope de estados a visualizar ($M$):** $(@bind M_slider Slider(4:1:12, default=8, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000063
md"**Estado inicial $Z_0$ para la cadena:** $(@bind z0_markov_slider Slider(1:1:8, default=2, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000064
begin
	M_dim = M_slider
	lambda_p = lambda_col_slider
	P_mat = zeros(Float64, M_dim + 1, M_dim + 1)
	
	# Fila 0: Estado absorbente
	P_mat[1, 1] = 1.0
	
	# Filas 1 a M: Construidas llamando a ModeloPBGW de src/markovModels.jl
	for i in 1:M_dim
		# Para i individuos con tasa individual lambda_p, la descendencia total sigue ModeloPBGW(i * lambda_p)
		m_i = ModeloPBGW(i * lambda_p; max_hijos = M_dim)
		for j in 0:(M_dim - 1)
			P_mat[i + 1, j + 1] = m_i.pk[j + 1]
		end
		# Todo el remanente de probabilidad va al estado agregado M+
		P_mat[i + 1, M_dim + 1] = max(0.0, 1.0 - sum(P_mat[i + 1, 1:M_dim]))
	end
	
	etiquetas_estados = ["$i" for i in 0:M_dim]
	etiquetas_estados[end] = "$M_dim+"
	
	# Paleta cálida y clara idéntica al modelo de parejas para garantizar máximo contraste con texto negro
	cmap_mat_markov = cgrad([:white, "#fef08a", "#fdba74", "#f87171", "#ef4444"], [0.0, 0.05, 0.35, 0.7, 1.0])
	fsize_p = M_dim <= 8 ? 10 : 8
	font_f = Sys.iswindows() ? "arialbd" : "sans-serif"
	
	p_heatmap_P = heatmap(
		etiquetas_estados, etiquetas_estados, P_mat,
		yflip = true,
		title = "Matriz de Transición P — Tamaño $(M_dim+1)x$(M_dim+1) (λ = $(@sprintf("%.2f", lambda_p)))",
		xlabel = "Estado destino en la siguiente generación (j)",
		ylabel = "Estado origen actual (i)",
		color = cmap_mat_markov, clims = (0.0, 1.0),
		size = (620, 520), aspect_ratio = :equal
	)
	
	for i in 1:(M_dim + 1)
		for j in 1:(M_dim + 1)
			v = P_mat[i, j]
			txt = @sprintf("%.2f", v)
			annotate!(p_heatmap_P, j, i, text(txt, font(fsize_p, :center, :black, font_f)))
		end
	end
	
	p_heatmap_P
end

# ╔═╡ 00000002-0002-0002-0002-000000000065
md"""
### Interpretación Directa de la Matriz:
- **Fila 0:** Es $[1.00, 0.00, 0.00, \ldots]$. Representa el **pozo absorbente**: una vez que el apellido llega a 0 portadores, nunca más sale.
- **Filas $i > 0$:** La probabilidad de pasar de $i$ personas a $j$ personas en la siguiente generación. 
- **Sesgo hacia la izquierda:** Como $\lambda = 0.5 < 1.0$, los valores más altos se acumulan a la izquierda de la diagonal (hacia números menores), empujando sistemáticamente la distribución hacia el estado 0.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000072
md"""
### Diagrama de Estados de la Cadena (Grafo Visual de Transición)

En un editor o simulador de Cadenas de Markov (como se diseña en teoría de autómatas o teoría de grafos), el proceso se visualiza como una **red dirigida**:
- **Nodos:** Representan los estados de la cadena (cantidad de portadores vivos). El estado **0** (rojo) es el pozo absorbente.
- **Flechas dirigidas:** Representan los saltos de transición en un paso (una generación), con su probabilidad $P_{ij}$.
- **Flechas inferiores rojas:** El camino de **decaimiento hacia la extinción** ($i \to j$ con $j < i$). Observa cómo son las más gruesas en régimen subcrítico ($\lambda = 0.5$).
- **Flechas superiores azules:** El camino de **crecimiento del linaje** ($i \to j$ con $j > i$).
- **Autolazo en 0:** Tiene probabilidad $1.00$ ($P_{0,0} = 1$). Si caes aquí, no hay flecha de salida.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000073
begin
	# Modelo de juguete para visualización de grafo: Estados 0, 1, 2, 3+
	# Construido directamente con ModeloPBGW de src/markovModels.jl
	M_toy = 3
	P_toy = zeros(Float64, M_toy + 1, M_toy + 1)
	P_toy[1, 1] = 1.0
	for i in 1:M_toy
		m_toy_i = ModeloPBGW(i * lambda_col_slider; max_hijos = M_toy)
		for j in 0:(M_toy - 1)
			P_toy[i + 1, j + 1] = m_toy_i.pk[j + 1]
		end
		P_toy[i + 1, M_toy + 1] = max(0.0, 1.0 - sum(P_toy[i + 1, 1:M_toy]))
	end
	
	xs_toy = [0.0, 3.2, 6.4, 9.6]
	ys_toy = [0.0, 0.0, 0.0, 0.0]
	nombres_toy = ["0\n(Extinto)", "1", "2", "3+"]
	colores_nodos = ["#e53935", "#1e88e5", "#1e88e5", "#8e24aa"]
	
	p_diag_markov = plot(
		size = (880, 430),
		xlims = (-1.8, 11.4), ylims = (-3.2, 3.4),
		aspect_ratio = :equal,
		title = "Grafo Dirigido de Transiciones de Estados (Ejemplo M=3, λ = $(@sprintf("%.2f", lambda_col_slider)))",
		legend = false, grid = false, framestyle = :none
	)
	
	r_nodo = 0.65
	
	# Función para dibujar arco con flecha orientada en dirección de transición
	function dibujar_flecha_curva!(plt, x1, y1, x2, y2, h_curva, col, lw, prob, es_superior)
		ctrl_x = (x1 + x2) / 2
		ctrl_y = 2 * h_curva
		
		# Punto t_arr donde la curva llega al borde exterior del nodo destino
		t_arr = 0.86
		for t_eval in range(0.70, 0.98, length=50)
			px = (1 - t_eval)^2 * x1 + 2 * (1 - t_eval) * t_eval * ctrl_x + t_eval^2 * x2
			py = (1 - t_eval)^2 * y1 + 2 * (1 - t_eval) * t_eval * ctrl_y + t_eval^2 * y2
			if hypot(px - x2, py - y2) <= r_nodo + 0.08
				t_arr = t_eval
				break
			end
		end
		
		# Tangente de la curva en t_arr para orientar la flecha
		tan_x = 2 * (1 - t_arr) * (ctrl_x - x1) + 2 * t_arr * (x2 - ctrl_x)
		tan_y = 2 * (1 - t_arr) * (ctrl_y - y1) + 2 * t_arr * (y2 - ctrl_y)
		L = hypot(tan_x, tan_y)
		ux = tan_x / L
		uy = tan_y / L
		
		tip_x = (1 - t_arr)^2 * x1 + 2 * (1 - t_arr) * t_arr * ctrl_x + t_arr^2 * x2
		tip_y = (1 - t_arr)^2 * y1 + 2 * (1 - t_arr) * t_arr * ctrl_y + t_arr^2 * y2
		
		arrow_len = 0.32
		arrow_w = 0.16
		base_x = tip_x - ux * arrow_len
		base_y = tip_y - uy * arrow_len
		
		ts_curva = range(0, t_arr, length=30)
		cx = (1 .- ts_curva).^2 .* x1 .+ 2 .* (1 .- ts_curva) .* ts_curva .* ctrl_x .+ ts_curva.^2 .* x2
		cy = (1 .- ts_curva).^2 .* y1 .+ 2 .* (1 .- ts_curva) .* ts_curva .* ctrl_y .+ ts_curva.^2 .* y2
		plot!(plt, cx, cy, color = col, lw = lw, alpha = 0.75)
		
		# Cabeza de flecha poligonal rellena apuntando exactamente al nodo destino
		left_x = base_x - uy * (arrow_w / 2)
		left_y = base_y + ux * (arrow_w / 2)
		right_x = base_x + uy * (arrow_w / 2)
		right_y = base_y - ux * (arrow_w / 2)
		plot!(plt, Shape([tip_x, left_x, right_x], [tip_y, left_y, right_y]),
		      color = col, linecolor = col, lw = 1, alpha = 0.95)
		
		mid_x = (x1 + x2) / 2
		mid_y = h_curva + (es_superior ? 0.28 : -0.28)
		annotate!(plt, mid_x, mid_y, text(@sprintf("%.2f", prob), 8, :center, col, Sys.iswindows() ? "arialbd" : "sans-serif"))
	end
	
	# Función para dibujar autolazo con flecha hacia el nodo
	function dibujar_autolazo_flecha!(plt, x0, y0, col, lw, prob)
		cx = x0
		cy = y0 + 1.05
		r_loop = 0.58
		thetas = range(-0.3, 2π + 0.3, length=50)
		lx = cx .+ r_loop .* cos.(thetas)
		ly = cy .+ r_loop .* sin.(thetas)
		plot!(plt, lx, ly, color = col, lw = lw, alpha = 0.75)
		
		# Punta de flecha en el lado izquierdo del bucle apuntando hacia abajo al nodo
		th = π + 0.5
		tip_x = cx + r_loop * cos(th)
		tip_y = cy + r_loop * sin(th)
		ux = sin(th)
		uy = -cos(th)
		L = hypot(ux, uy)
		ux /= L; uy /= L
		
		al = 0.28; aw = 0.16
		bx = tip_x - ux * al
		by = tip_y - uy * al
		plot!(plt, Shape([tip_x, bx - uy * (aw/2), bx + uy * (aw/2)], [tip_y, by + ux * (aw/2), by - ux * (aw/2)]),
		      color = col, linecolor = col, lw = 1, alpha = 0.95)
		
		annotate!(plt, cx, cy + r_loop + 0.32, text(@sprintf("%.2f", prob), 8, :center, col, Sys.iswindows() ? "arialbd" : "sans-serif"))
	end
	
	# Dibujar todas las transiciones con flechas según P_toy
	for i in 0:M_toy
		for j in 0:M_toy
			prob = P_toy[i + 1, j + 1]
			if prob > 0.025
				x_i, y_i = xs_toy[i + 1], ys_toy[i + 1]
				x_j, y_j = xs_toy[j + 1], ys_toy[j + 1]
				
				if i == j
					col_l = i == 0 ? "#e53935" : "#757575"
					lw_l = i == 0 ? 2.5 : 1.5
					dibujar_autolazo_flecha!(p_diag_markov, x_i, y_i, col_l, lw_l, prob)
				elseif j > i
					# Transición creciente i -> j (hacia la derecha, arco superior azul)
					dist_step = abs(j - i)
					h_val = 0.55 + 0.40 * (dist_step - 1)
					col_u = "#1e88e5"
					lw_u = 1.2 + 3.0 * prob
					dibujar_flecha_curva!(p_diag_markov, x_i, y_i, x_j, y_j, h_val, col_u, lw_u, prob, true)
				else
					# Transición decreciente i -> j (hacia la izquierda / extinción, arco inferior rojo/naranja)
					dist_step = abs(i - j)
					h_val = -(0.65 + 0.55 * (dist_step - 1))
					col_d = j == 0 ? "#e53935" : "#e65100"
					lw_d = 1.2 + 3.5 * prob
					dibujar_flecha_curva!(p_diag_markov, x_i, y_i, x_j, y_j, h_val, col_d, lw_d, prob, false)
				end
			end
		end
	end
	
	# Dibujar nodos circulares con etiquetas
	for idx in 1:(M_toy + 1)
		scatter!(p_diag_markov, [xs_toy[idx]], [ys_toy[idx]],
		         markersize = 28, color = :white,
		         markerstrokecolor = colores_nodos[idx],
		         markerstrokewidth = 3.5)
		annotate!(p_diag_markov, xs_toy[idx], ys_toy[idx],
		          text(nombres_toy[idx], 9, :center, :black, Sys.iswindows() ? "arialbd" : "sans-serif"))
	end
	
	p_diag_markov
end

# ╔═╡ 00000002-0002-0002-0002-000000000066
md"**Paso de tiempo / Generación ($n$) para el vector de estado $\boldsymbol{\pi}_n$:** $(@bind gen_pi_slider Slider(0:1:25, default=3, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000067
begin
	M_curr = M_slider
	z0_c = min(z0_markov_slider, M_curr)
	
	# Vector inicial: 100% en z0_c
	pi_0 = zeros(Float64, M_curr + 1)
	pi_0[z0_c + 1] = 1.0
	
	# Calcular pi_n = pi_0 * P^n
	pi_n = copy(pi_0)
	for _ in 1:gen_pi_slider
		pi_n = pi_n' * P_mat
		pi_n = vec(pi_n)
	end
	
	colores_pi = [k == 0 ? "#e53935" : "#1e88e5" for k in 0:M_curr]
	
	p_bar_pi = bar(
		etiquetas_estados, pi_n,
		xlabel = "Número de Portadores (Estado)",
		ylabel = "Probabilidad πₙ(k)",
		title = "Distribución de Probabilidad en la Generación n = $(gen_pi_slider) (Arrancando con Z₀ = $(z0_c))",
		color = colores_pi, alpha = 0.85, legend = false,
		ylims = (0.0, 1.05), size = (800, 380)
	)
	
	for (idx, p_val) in enumerate(pi_n)
		if p_val > 0.015
			annotate!(p_bar_pi, idx, p_val + 0.03, text(@sprintf("%.1f%%", p_val * 100), 8, :center, :black, Sys.iswindows() ? "arialbd" : "sans-serif"))
		end
	end
	
	p_bar_pi
end

# ╔═╡ 00000002-0002-0002-0002-000000000068
begin
	max_g_flujo = 20
	M_flujo = M_slider
	z0_f = min(z0_markov_slider, M_flujo)
	
	mat_flujo = zeros(Float64, max_g_flujo + 1, M_flujo + 1)
	vec_temp = zeros(Float64, M_flujo + 1)
	vec_temp[z0_f + 1] = 1.0
	mat_flujo[1, :] = vec_temp
	
	for g in 1:max_g_flujo
		vec_temp = vec(vec_temp' * P_mat)
		mat_flujo[g + 1, :] = vec_temp
	end
	
	etiquetas_gens = ["Gen $g" for g in 0:max_g_flujo]
	
	p_flujo = heatmap(
		etiquetas_estados, etiquetas_gens, mat_flujo,
		yflip = true,
		title = "Evolución Temporal del Flujo de Probabilidad: De Z₀ = $(z0_f) a la Extinción",
		xlabel = "Estado (Portadores)", ylabel = "Generación (Paso temporal)",
		color = cmap_mat_markov, clims = (0.0, 1.0),
		size = (700, 450)
	)
	p_flujo
end

# ╔═╡ 00000002-0002-0002-0002-000000000069
begin
	gens_curva_p = 0:25
	prob_ext_teorica = Float64[]
	vec_t = zeros(Float64, M_dim + 1)
	vec_t[min(z0_markov_slider, M_dim) + 1] = 1.0
	
	for g in gens_curva_p
		push!(prob_ext_teorica, vec_t[1])
		vec_t = vec(vec_t' * P_mat)
	end
	
	p_ext_acum = plot(
		gens_curva_p, prob_ext_teorica,
		xlabel = "Generación (n)", ylabel = "Probabilidad acumulada P(Zₙ = 0)",
		title = "Curva de Absorción: ¿Con qué rapidez la cadena cae al estado 0?",
		label = "P(Extinción acumulada)", color = :red, linewidth = 3,
		size = (800, 360), ylims = (0.0, 1.05), grid = true
	)
	hline!(p_ext_acum, [0.5], label = "Vida media (50%)", color = :gray, linestyle = :dash)
	hline!(p_ext_acum, [1.0], label = "Absorción total (100%)", color = :black, linestyle = :dot)
	p_ext_acum
end

# ╔═╡ 00000002-0002-0002-0002-000000000070
md"""
### Caracterización Matemática Formal de la Cadena

| Propiedad | Clasificación en este Modelo | Explicación Simplificada |
|:---|:---|:---|
| **Espacio de Estados** | $\{0, 1, 2, \ldots, M\}$ | Número de personas vivas portadoras del apellido |
| **Clasificación del Estado 0** | **Absorbente y Recurrente** | Si se extingue, $P_{0,0} = 1$. No vuelve a surgir |
| **Clasificación de Estados $k \ge 1$** | **Transitorios** | La probabilidad de volver a ellos tiende a 0 conforme $n \to \infty$ |
| **Irreducibilidad** | **No irreducible** | No todos los estados se comunican: de los transitorios vas al 0, pero del 0 no puedes volver |
| **Periodicidad** | **Aperiódica** ($d = 1$) | El estado 0 tiene bucle $P_{0,0} = 1$; no hay ciclos oscilatorios forzados |
| **Distribución Límite / Estacionaria** | $\boldsymbol{\pi}^* = [1.0, 0, 0, \ldots, 0]$ | Toda la masa probabilística termina inevitablemente concentrada en el estado 0 |
| **Tiempo de Absorción** | Finito casi seguramente | En régimen subcrítico ($\lambda < 1$), el tiempo hasta la extinción es finito con probabilidad 1 |
"""

# ╔═╡ 00000002-0002-0002-0002-000000000071
md"""
### Resumen simplificado

Imagina la cadena de Markov de Galton-Watson como un modelo dinámico paso a paso:
1. **El espacio de estados:** Cada posición representa un número de portadores ($0, 1, 2, 3 \ldots$).
2. **El estado inicial:** Arranca en la cantidad actual del apellido ($Z_0$).
3. **La matriz de transición $P$:** En cada generación, la población salta a un nuevo estado según las probabilidades calculadas en su fila actual.
4. **El régimen subcrítico ($\lambda < 1$):** Las probabilidades se concentran a la izquierda de la diagonal, generando una deriva sistemática hacia números menores.
5. **El estado 0 (Extinción):** Es un pozo absorbente. Una vez que se alcanza 0, no hay descendencia posible y la cadena permanece en 0 con probabilidad 1.
6. **Conclusión estocástica:** Sin importar si el apellido cuenta con 1 individuo o con los 754.991 Rodríguez del Censo DANE, con fertilidad por debajo del reemplazo, **toda la probabilidad converge inexorablemente hacia el estado 0**.
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
# ╟─00000002-0002-0002-0002-000000000014
# ╟─00000002-0002-0002-0002-000000000015
# ╟─00000002-0002-0002-0002-000000000016
# ╟─00000002-0002-0002-0002-000000000017
# ╟─00000002-0002-0002-0002-000000000018
# ╟─00000002-0002-0002-0002-000000000019
# ╟─00000002-0002-0002-0002-000000000020
# ╟─00000002-0002-0002-0002-000000000021
# ╟─00000002-0002-0002-0002-000000000022
# ╟─00000002-0002-0002-0002-000000000023
# ╟─00000002-0002-0002-0002-000000000024
# ╟─00000002-0002-0002-0002-000000000025
# ╟─00000002-0002-0002-0002-000000000026
# ╟─00000002-0002-0002-0002-000000000027
# ╟─00000002-0002-0002-0002-000000000028
# ╟─00000002-0002-0002-0002-000000000029
# ╟─00000002-0002-0002-0002-000000000030
# ╟─00000002-0002-0002-0002-000000000031
# ╟─00000002-0002-0002-0002-000000000032
# ╟─00000002-0002-0002-0002-000000000040
# ╟─00000002-0002-0002-0002-000000000041
# ╟─00000002-0002-0002-0002-000000000042
# ╟─00000002-0002-0002-0002-000000000043
# ╟─00000002-0002-0002-0002-000000000044
# ╟─00000002-0002-0002-0002-000000000045
# ╟─00000002-0002-0002-0002-000000000046
# ╟─00000002-0002-0002-0002-000000000047
# ╟─00000002-0002-0002-0002-000000000048
# ╟─00000002-0002-0002-0002-000000000049
# ╟─00000002-0002-0002-0002-000000000050
# ╟─00000002-0002-0002-0002-000000000051
# ╟─00000002-0002-0002-0002-000000000052
# ╟─00000002-0002-0002-0002-000000000060
# ╟─00000002-0002-0002-0002-000000000061
# ╟─00000002-0002-0002-0002-000000000062
# ╟─00000002-0002-0002-0002-000000000063
# ╟─00000002-0002-0002-0002-000000000064
# ╟─00000002-0002-0002-0002-000000000065
# ╟─00000002-0002-0002-0002-000000000072
# ╟─00000002-0002-0002-0002-000000000073
# ╟─00000002-0002-0002-0002-000000000066
# ╟─00000002-0002-0002-0002-000000000067
# ╟─00000002-0002-0002-0002-000000000068
# ╟─00000002-0002-0002-0002-000000000069
# ╟─00000002-0002-0002-0002-000000000070
# ╟─00000002-0002-0002-0002-000000000071
