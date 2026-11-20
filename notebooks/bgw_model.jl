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
md"**Simulaciones Monte Carlo:** $(@bind n_sims_col_slider Slider(100:100:2000, default=500, show_value=true))"

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
	
	# Tabla comparativa resaltando el apellido activo
	filas_tabla = [
		"| Puesto | Apellido | Portadores DANE 2018 ($Z_0$) | Estado en Simulación |",
		"|:---|:---|:---|:---|"
	]
	for (idx, (ap, cnt)) in enumerate(zip(apellidos_colombia, conteos_dane))
		marca = ap == ap_col_select ? "**👉 SELECCIONADO**" : "En espera"
		ap_fmt = ap == ap_col_select ? "**$ap**" : ap
		push!(filas_tabla, "| $idx | $ap_fmt | $(cnt) | $marca |")
	end
	Markdown.parse(join(filas_tabla, "\n"))
end

# ╔═╡ 00000002-0002-0002-0002-000000000047
md"""
### ⚙️ Trucos Numéricos y Tipo de Monte Carlo Utilizado

#### 1. El Truco Numérico Aditivo de Poisson ($O(1)$ vs $O(N)$)
En una implementación ingenua de Galton-Watson, para una generación con $Z_n = 754.991$ individuos, tendríamos que generar 754.991 números aleatorios en cada generación. Con 30 generaciones y 500 simulaciones, ¡eso requeriría generar más de **11 mil millones de variables aleatorias**, desbordando la memoria RAM y tardando horas!

**El truco matemático en Julia:**
Aprovechamos la **propiedad reproductiva de la distribución de Poisson**: si $X_1, X_2, \ldots, X_{Z_n} \stackrel{iid}{\sim} \text{Poisson}(\lambda)$, entonces su suma exacta es:
$$\sum_{i=1}^{Z_n} X_i \sim \text{Poisson}(Z_n \cdot \lambda)$$
Por lo tanto, en Julia podemos simular la generación completa en **una sola llamada $O(1)$**:
```julia
z = rand(rng, Poisson(z * lambda))
```
Esto reduce el tiempo de cómputo de horas a **milisegundos**, permitiendo correr 2.000 simulaciones completas en tiempo real.

#### 2. ¿Qué tipo de Monte Carlo usamos? (Forward MC vs MCMC)
- **Monte Carlo Directo / Hacia Adelante (Forward Monte Carlo):** Es el que utilizamos aquí. Muestreamos realizaciones estocásticas directas del proceso temporal $\{Z_n\}_{n \ge 0}$ a lo largo de las generaciones.
- **¿Por qué NO usamos MCMC (Markov Chain Monte Carlo - Metropolis/Gibbs)?**
  MCMC se usa cuando queremos muestrear de una distribución de probabilidad estacionaria o a posteriori muy difícil de calcular de forma analítica, creando una cadena artificial que converge a ella. En Galton-Watson subcrítico, el proceso es **transitorio y absorbente** (termina en 0), por lo que lo que nos interesa es la **dinámica temporal de decaimiento**, la cual se modela naturalmente con simulación directa hacia adelante (Forward MC).
"""

# ╔═╡ 00000002-0002-0002-0002-000000000048
begin
	# Función optimizada con el truco numérico aditivo
	function simular_extincion_colombia(lambda, z0, max_gen, rng)
		z = z0
		gen_ext = max_gen + 1
		trayectoria = Int[z0]
		
		for g in 1:max_gen
			if z == 0
				push!(trayectoria, 0)
				if gen_ext > max_gen
					gen_ext = g - 1
				end
				continue
			end
			# Propiedad aditiva de Poisson O(1)
			z = rand(rng, Poisson(z * lambda))
			push!(trayectoria, z)
			if z == 0 && gen_ext > max_gen
				gen_ext = g
			end
		end
		return (extinto = (trayectoria[end] == 0), gen_ext = gen_ext, trayectoria = trayectoria)
	end

	rng_mc_col = MersenneTwister(semilla_col_slider)
	
	max_gen_sim = max(gen_eval_slider + 15, 45)
	tiempos_ext = Int[]
	trayectorias_muestra = Vector{Vector{Int}}()
	extintos_en_X = 0
	
	for s in 1:n_sims_col_slider
		res = simular_extincion_colombia(lambda_col_slider, z0_seleccionado, max_gen_sim, rng_mc_col)
		if s <= 6
			push!(trayectorias_muestra, res.trayectoria)
		end
		if res.gen_ext <= max_gen_sim
			push!(tiempos_ext, res.gen_ext)
		end
		if res.gen_ext <= gen_eval_slider
			extintos_en_X += 1
		end
	end
	
	prob_ext_X = extintos_en_X / n_sims_col_slider
	mediana_ext = isempty(tiempos_ext) ? 0 : round(Int, median(tiempos_ext))
	media_ext = isempty(tiempos_ext) ? 0.0 : mean(tiempos_ext)
	
	# Estimación teórica continua aproximada: T ≈ ln(Z0) / ln(1/lambda)
	t_teorico = lambda_col_slider < 1.0 ? log(z0_seleccionado) / log(1.0 / lambda_col_slider) : Inf
	
	md"""
	### 📊 Resultados del Monte Carlo para **$(ap_col_select)**
	
	| Métrica | Valor Estimado | Significado |
	|:---|:---|:---|
	| **Población Inicial ($Z_0$)** | **$(z0_seleccionado)** | Portadores según Censo DANE 2018 |
	| **Tasa de natalidad masculina ($\lambda$)** | **$(lambda_col_slider)** | $(lambda_col_slider < 1.0 ? "🔴 Subcrítico (Extinción inevitable)" : (lambda_col_slider == 1.0 ? "🟡 Crítico" : "🟢 Supercrítico")) |
	| **Probabilidad de extinción en $\le $(gen_eval_slider) generaciones** | <span style="font-size:1.3em; font-weight:bold; color:$(prob_ext_X > 0.5 ? "crimson" : "steelblue");">$(round(prob_ext_X * 100, digits=2))%</span> | Fracción de simulaciones extintas en Gen $(gen_eval_slider) |
	| **Generación Mediana de Extinción** | **Gen $(mediana_ext)** (~$(mediana_ext * 28) años) | El 50% de las historias se extinguieron antes de este año |
	| **Generación Media de Extinción** | **Gen $(round(media_ext, digits=1))** | Promedio aritmético de las simulaciones |
	| **Aproximación Teórica Determinista** | **Gen $(round(t_teorico, digits=1))** | Fórmula: $T \approx \frac{\ln(Z_0)}{\ln(1/\lambda)}$ |
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
	# Comparativa de los 10 apellidos
	tiempos_medios_10 = Float64[]
	rng_comp = MersenneTwister(123)
	for cnt in conteos_dane
		tiempos_i = Int[]
		for _ in 1:250
			res_i = simular_extincion_colombia(lambda_col_slider, cnt, 45, rng_comp)
			push!(tiempos_i, res_i.gen_ext)
		end
		push!(tiempos_medios_10, median(tiempos_i))
	end
	
	colores_barras = [ap == ap_col_select ? "#e53935" : "#1e88e5" for ap in apellidos_colombia]
	
	p_comp = bar(
		apellidos_colombia, tiempos_medios_10,
		orientation = :horizontal,
		xlabel = "Generaciones estimadas hasta extinción",
		ylabel = "Apellido",
		title = "¿Cuánto tarda en extinguirse cada apellido? (λ = $(lambda_col_slider))",
		color = colores_barras, alpha = 0.85, legend = false,
		size = (800, 420), xlims = (0, maximum(tiempos_medios_10) * 1.25)
	)
	
	for (idx, (ap, t_med)) in enumerate(zip(apellidos_colombia, tiempos_medios_10))
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
### 🏷️ 1. Definición del Espacio de Estados
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
	
	# Filas 1 a M: Estados transitorios
	for i in 1:M_dim
		pois_i = Poisson(i * lambda_p)
		for j in 0:(M_dim - 1)
			P_mat[i + 1, j + 1] = pdf(pois_i, j)
		end
		# Todo el remanente de probabilidad va al estado agregado M
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
### 🔍 Interpretación Directa de la Matriz:
- **Fila 0:** Es $[1.00, 0.00, 0.00, \ldots]$. Representa el **pozo absorbente**: una vez que el apellido llega a 0 portadores, nunca más sale.
- **Filas $i > 0$:** La probabilidad de pasar de $i$ personas a $j$ personas en la siguiente generación. 
- **Sesgo hacia la izquierda:** Como $\lambda = 0.5 < 1.0$, los valores más altos se acumulan a la izquierda de la diagonal (hacia números menores), empujando sistemáticamente la distribución hacia el estado 0.
"""

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
### 📋 Caracterización Matemática Formal de la Cadena

| Propiedad | Clasificación en este Modelo | Explicación en Cristiano |
|:---|:---|:---|
| **Espacio de Estados** | $\{0, 1, 2, \ldots, M\}$ | Número de personas vivas portadoras del apellido |
| **Clasificación del Estado 0** | **Absorbente y Recurrente** | Si se extingue, $P_{0,0} = 1$. Nunca resucita |
| **Clasificación de Estados $k \ge 1$** | **Transitorios** | La probabilidad de volver a ellos tiende a 0 conforme $n \to \infty$ |
| **Irreducibilidad** | ❌ **No irreducible** | No todos los estados se comunican: de los transitorios vas al 0, pero del 0 no puedes volver |
| **Periodicidad** | **Aperiódica** ($d = 1$) | El estado 0 tiene bucle $P_{0,0} = 1$; no hay ciclos oscilatorios forzados |
| **Distribución Límite / Estacionaria** | $\boldsymbol{\pi}^* = [1.0, 0, 0, \ldots, 0]$ | Toda la masa probabilística termina inevitablemente concentrada en el estado 0 |
| **Tiempo de Absorción** | Finito casi seguramente | En régimen subcrítico ($\lambda < 1$), el tiempo hasta la extinción es finito con probabilidad 1 |
"""

# ╔═╡ 00000002-0002-0002-0002-000000000071
md"""
### 💡 Resumen para Humanos ("Para Dummies")

Imagina la cadena de Markov de Galton-Watson como un **juego de mesa con un tobogán empinado**:
1. **El tablero:** Cada casilla es un número de personas ($0, 1, 2, 3 \ldots$).
2. **Tu ficha:** Empieza en la casilla de tu apellido ($Z_0$).
3. **El dado (la matriz $P$):** En cada turno (generación), tiras los dados según la fila de tu casilla actual.
4. **La trampa del régimen subcrítico ($\lambda < 1$):** El dado está cargado para dar números menores. El tablero está inclinado hacia abajo.
5. **La casilla 0:** Es un agujero negro. Una vez que caes en 0, no hay dados que tirar: te quedas en 0 para siempre.
6. **El resultado:** No importa si empiezas con 1 persona o con los 754.991 Rodríguez de Colombia; si la fertilidad es menor a la de reemplazo, **toda la probabilidad resbala por el tobogán y termina dentro del hoyo 0**.
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
# ╟─00000002-0002-0002-0002-000000000066
# ╟─00000002-0002-0002-0002-000000000067
# ╟─00000002-0002-0002-0002-000000000068
# ╟─00000002-0002-0002-0002-000000000069
# ╟─00000002-0002-0002-0002-000000000070
# ╟─00000002-0002-0002-0002-000000000071
