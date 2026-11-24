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

# ╔═╡ d0000003-0003-0003-0003-000000000001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", ".."))
    using Random, Distributions, Plots, PlutoUI, Printf, StatsBase

    # ── Modelo BGW y datos reales, cargados tal cual de src/ ──
    include(joinpath(@__DIR__, "..", "..", "src", "core", "numericalMethods.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "bgwModel.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "experiments", "daneData.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "experiments", "daneBgwRunners.jl"))
    md"Entorno, modelo BGW y datos DANE cargados desde `src/`."
end

# ╔═╡ d0000003-0003-0003-0003-000000000002
md"""
# Resultados: ¿Cuándo se extinguirán los apellidos de Colombia?
## Proceso de Bienaymé-Galton-Watson con datos del Censo DANE 2018

Este notebook aplica el modelo BGW (`src/bgwModel.jl`) a los 10 apellidos más
frecuentes de Colombia (Censo Nacional de Población y Vivienda, DANE 2018). Los
datos y la orquestación Monte Carlo viven en `src/experiments/daneData.jl`; aquí
solo elegimos parámetros y dibujamos.

### Contexto demográfico

La tasa de fecundidad en Colombia ha caído por debajo del reemplazo. Como el
primer apellido tradicionalmente lo transmiten los hijos varones, cada hombre
engendra en promedio $\lambda \approx 0.5$ hijos varones que continúan el
linaje. Con $m = \lambda \le 1$ el proceso es **subcrítico**: la extinción es
matemáticamente segura ($q = 1$). La pregunta deja de ser *si* se extingue y
pasa a ser *en cuántas generaciones*.
"""

# ╔═╡ d0000003-0003-0003-0003-000000000010
@bind idx_apellido Select([i => "$i. $(APELLIDOS_DANE[i]) ($(CONTEOS_DANE[i]) personas)" for i in 1:10], default=1)

# ╔═╡ d0000003-0003-0003-0003-000000000011
md"**Tasa de natalidad masculina $\lambda$** (hijos varones promedio por hombre): $(@bind lambda_col Slider(0.1:0.05:1.2, default=0.5, show_value=true))"

# ╔═╡ d0000003-0003-0003-0003-000000000012
md"**Evaluar probabilidad de extinción acumulada a las generaciones ($X$):** $(@bind gen_eval Slider(5:1:40, default=20, show_value=true))"

# ╔═╡ d0000003-0003-0003-0003-000000000013
md"**Simulaciones Monte Carlo:** $(@bind n_sims_col Slider(50:25:500, default=100, show_value=true))"

# ╔═╡ d0000003-0003-0003-0003-000000000014
md"**Semilla aleatoria:** $(@bind semilla_col Slider(1:100, default=42, show_value=true))"

# ╔═╡ d0000003-0003-0003-0003-000000000020
begin
    ap_idx = idx_apellido
    ap_nombre = APELLIDOS_DANE[ap_idx]
    z0_sel = CONTEOS_DANE[ap_idx]

    # Modelo BGW con el constructor de src/bgwModel.jl.
    modelo_col = ModeloPBGW(lambda_col; max_hijos = 15)

    filas_tabla = [
        "| Puesto | Apellido | Portadores DANE 2018 (\$Z_0\$) | Estado |",
        "|:---|:---|:---|:---|",
    ]
    for i in 1:10
        marca = i == ap_idx ? "**SELECCIONADO**" : "En espera"
        ap_fmt = i == ap_idx ? "**$(APELLIDOS_DANE[i])**" : APELLIDOS_DANE[i]
        push!(filas_tabla, "| $i | $ap_fmt | $(CONTEOS_DANE[i]) | $marca |")
    end
    Markdown.parse(join(filas_tabla, "\n"))
end

# ╔═╡ d0000003-0003-0003-0003-000000000021
begin
    # Monte Carlo del tiempo de extincion: runner de src/experiments/daneData.jl.
    res_col = correr_montecarlo_colombia(modelo_col, z0_sel, n_sims_col, gen_eval, semilla_col)

    # Validacion teorica con funciones de src/bgwModel.jl.
    m_esp = esperanza(modelo_col)
    es_segura = es_extincion_segura(modelo_col)
    q_teorico = probabilidad_extincion(modelo_col)
    t_teorico = m_esp < 1.0 ? log(z0_sel) / log(1.0 / m_esp) : Inf

    md"""
    ### Resultados del Monte Carlo para **$(ap_nombre)**

    | Métrica | Valor | Función de `src/` |
    |:---|:---|:---|
    | Población inicial $Z_0$ | **$(z0_sel)** | Dato DANE 2018 |
    | Esperanza $m = E[X]$ | **$(round(m_esp, digits=3))** | `esperanza` |
    | Régimen | **$(es_segura ? "Subcrítico/Crítico (m ≤ 1)" : "Supercrítico (m > 1)")** | `es_extincion_segura` |
    | Probabilidad asintótica $q$ | **$(round(q_teorico * 100, digits=2))%** | `probabilidad_extincion` |
    | P(extinción en ≤ $(gen_eval) gen) | **$(round(res_col.prob_ext_X * 100, digits=2))%** | `correr_montecarlo_colombia` |
    | Generación mediana de extinción | **Gen $(res_col.mediana_ext)** (~$(res_col.mediana_ext * 28) años) | mediana Monte Carlo |
    | Generación media de extinción | **Gen $(round(res_col.media_ext, digits=1))** | media Monte Carlo |
    | Aproximación determinista | **Gen $(round(t_teorico, digits=1))** | $T \approx \ln(Z_0)/\ln(1/m)$ |
    """
end

# ╔═╡ d0000003-0003-0003-0003-000000000022
begin
    # Histograma del tiempo de extincion.
    if !isempty(res_col.tiempos_ext)
        p_hist = histogram(
            res_col.tiempos_ext,
            bins = range(minimum(res_col.tiempos_ext) - 0.5, maximum(res_col.tiempos_ext) + 0.5, step = 1),
            xlabel = "Generación de extinción (1 gen ≈ 25-30 años)",
            ylabel = "Frecuencia de simulaciones",
            title = "Tiempo de extinción para $(ap_nombre) (n=$(n_sims_col))",
            color = "#e53935", alpha = 0.75, legend = :topright, size = (800, 380),
        )
        vline!(p_hist, [res_col.mediana_ext], label = "Mediana = Gen $(res_col.mediana_ext)", color = :black, lw = 2.5, linestyle = :dash)
        vline!(p_hist, [gen_eval], label = "Umbral X = Gen $(gen_eval) ($(round(res_col.prob_ext_X*100, digits=1))% extinto)", color = :blue, lw = 2)
        p_hist
    else
        plot(title = "No hubo extinciones (λ ≥ 1.0)")
    end
end

# ╔═╡ d0000003-0003-0003-0003-000000000023
begin
    # Trayectorias en escala logaritmica.
    p_log = plot(
        xlabel = "Generación", ylabel = "Portadores (escala log10)",
        title = "Decaimiento estocástico — $(ap_nombre) (λ = $(@sprintf("%.2f", lambda_col)))",
        yaxis = :log10, size = (800, 400), legend = :bottomleft,
        grid = true, gridstyle = :dash, gridalpha = 0.3,
    )
    gens_curva = 0:min(res_col.max_gen_sim, 35)
    curva_esp = [max(1.0, z0_sel * (lambda_col^g)) for g in gens_curva]
    plot!(p_log, gens_curva, curva_esp, label = "Media esperada E[Zₙ] = Z₀·λⁿ", color = :black, lw = 2.5, linestyle = :dash)
    for (idx, tr) in enumerate(res_col.trayectorias_muestra)
        tr_plot = [max(0.5, Float64(x)) for x in tr]
        plot!(p_log, 0:(length(tr_plot) - 1), tr_plot, label = idx == 1 ? "Simulaciones individuales" : false, alpha = 0.7, lw = 1.5)
    end
    hline!(p_log, [1.0], label = "1 portador (umbral crítico)", color = :red, linestyle = :dot, lw = 1.5)
    p_log
end

# ╔═╡ d0000003-0003-0003-0003-000000000024
md"""
### Comparativa de los 10 apellidos

Generación mediana de extinción de cada apellido (poblaciones iniciales del
Censo DANE 2018), estimada con `calcular_medianas_10` de `src/experiments/daneData.jl`.
"""

# ╔═╡ d0000003-0003-0003-0003-000000000025
begin
    tiempos_medios_10 = calcular_medianas_10(modelo_col)
    colores_barras = [i == ap_idx ? "#e53935" : "#1e88e5" for i in 1:10]

    aps_rev = reverse(APELLIDOS_DANE)
    tiempos_rev = reverse(tiempos_medios_10)
    colores_rev = reverse(colores_barras)

    p_comp = bar(
        1:10, tiempos_rev,
        orientation = :horizontal,
        yticks = (1:10, aps_rev),
        xlabel = "Generaciones estimadas hasta extinción",
        ylabel = "Apellido",
        title = "¿Cuánto tarda en extinguirse cada apellido? (λ = $(@sprintf("%.2f", lambda_col)))",
        color = colores_rev, alpha = 0.85, legend = false,
        size = (800, 420), xlims = (0, maximum(tiempos_medios_10) * 1.35),
    )
    for (idx, t_med) in enumerate(tiempos_rev)
        annotate!(p_comp, t_med + 0.4, idx, text("Gen $(round(Int, t_med)) (~$(round(Int, t_med*28)) años)", 8, :left, :black))
    end
    p_comp
end

# ╔═╡ d0000003-0003-0003-0003-000000000026
md"""
> **La paradoja de los apellidos masivos:** en un proceso subcrítico la población
> decae geométricamente ($Z_n \approx Z_0 \lambda^n$), así que el tiempo de
> extinción crece con el **logaritmo** de la población inicial:
> $$T_{\text{ext}} \approx \frac{\ln(Z_0)}{\ln(1/\lambda)}$$
> Duplicar el número de portadores apenas añade **una generación** de vida al
> apellido. Por eso Rodríguez (754.991) y Ramírez (441.179) se extinguen casi a
> la vez pese a los 300.000 portadores de diferencia.
"""

# ╔═╡ Cell order:
# ╟─d0000003-0003-0003-0003-000000000001
# ╟─d0000003-0003-0003-0003-000000000002
# ╠═d0000003-0003-0003-0003-000000000010
# ╟─d0000003-0003-0003-0003-000000000011
# ╟─d0000003-0003-0003-0003-000000000012
# ╟─d0000003-0003-0003-0003-000000000013
# ╟─d0000003-0003-0003-0003-000000000014
# ╟─d0000003-0003-0003-0003-000000000020
# ╟─d0000003-0003-0003-0003-000000000021
# ╠═d0000003-0003-0003-0003-000000000022
# ╠═d0000003-0003-0003-0003-000000000023
# ╟─d0000003-0003-0003-0003-000000000024
# ╠═d0000003-0003-0003-0003-000000000025
# ╟─d0000003-0003-0003-0003-000000000026
