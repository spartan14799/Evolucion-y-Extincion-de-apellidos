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

# ╔═╡ e0000004-0004-0004-0004-000000000001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", ".."))
    using Random, StatsBase, Plots, PlutoUI, Printf

    # ── Modelo por parejas y datos reales, cargados tal cual de src/ ──
    include(joinpath(@__DIR__, "..", "..", "src", "coupleModel.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "core", "graphTopologies.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "experiments", "daneData.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "experiments", "daneCoupleRunners.jl"))
    md"Entorno, modelo por parejas y datos DANE cargados desde `src/`."
end

# ╔═╡ e0000004-0004-0004-0004-000000000002
md"""
# Resultados: árbol genealógico de apellidos (Censo DANE 2018)

Simulamos el árbol genealógico generación por generación con el motor del modelo
por parejas (`siguiente_generacion_parejas` de `src/coupleModel.jl`) sobre los 10
apellidos más frecuentes de Colombia (DANE 2018). La simulación la orquesta
`simular_arbol_dane` de `src/experiments/daneData.jl`; este notebook solo dibuja.

Cada fila es una generación, cada círculo una persona (el color es su apellido)
y las líneas conectan cada pareja de padres con sus dos hijos. La unidad que se
reproduce es la **pareja**: cada pareja de la generación `t` aporta dos hijos a
la generación `t+1`, que heredan el apellido de uno de los dos padres.
"""

# ╔═╡ e0000004-0004-0004-0004-000000000010
md"""**Modo de población inicial**: $(@bind modo_dane Select([
    "equitativo" => "Equitativo (reparte parejas entre los 10 apellidos)",
    "dane" => "Proporcional al Censo DANE 2018",
], default="equitativo"))"""

# ╔═╡ e0000004-0004-0004-0004-000000000011
md"""**Parejas totales** (población = 2 × parejas): $(@bind parejas_dane Slider(5:1:20, default=10, show_value=true))"""

# ╔═╡ e0000004-0004-0004-0004-000000000012
md"""**Modo de probabilidad**: $(@bind modo_prob_dane Select([
    "ponderado" => "Ponderado por conteo",
    "uniforme_aristas" => "Uniforme por aristas",
], default="ponderado"))"""

# ╔═╡ e0000004-0004-0004-0004-000000000013
md"""**Semilla de la partida**: $(@bind semilla_dane Slider(1:100, default=7, show_value=true))"""

# ╔═╡ e0000004-0004-0004-0004-000000000014
md"""**Máximo de generaciones a mostrar**: $(@bind max_gen_dane Slider(5:5:40, default=25, show_value=true))"""

# ╔═╡ e0000004-0004-0004-0004-000000000020
begin
    # Simulacion del arbol con el runner de src/experiments/daneData.jl.
    arbol = simular_arbol_dane(;
        modo_inicial = Symbol(modo_dane),
        parejas = parejas_dane,
        modo_prob = Symbol(modo_prob_dane),
        semilla = semilla_dane,
        max_generaciones = max_gen_dane,
    )
    md"Simulación lista: $(arbol.generaciones) generaciones, $(arbol.n_total) personas."
end

# ╔═╡ e0000004-0004-0004-0004-000000000021
begin
    # Dibujo del arbol: una columna por persona, una fila por generacion.
    hist = arbol.historial
    G = arbol.generaciones
    N = arbol.n_total
    rng_dibujo = MersenneTwister(semilla_dane + 1000)

    p_arbol = plot(size = (820, max(460, 30 * (G + 1))),
                   legend = :outerright,
                   xlims = (0.3, N + 0.7),
                   ylims = (-G - 0.7, 0.7),
                   yticks = (-(0:G), ["Gen $g" for g in 0:G]),
                   title = "Árbol genealógico por parejas (DANE 2018)",
                   grid = false, framestyle = :box)

    # Lineas de parentesco: los hijos vienen agrupados de dos en dos por pareja;
    # conectamos cada par de hijos con un padre de su mismo apellido.
    xc = Float64[]; yc = Float64[]
    for g in 1:G
        padres = hist[g]
        hijos = hist[g + 1]
        np = length(hijos)
        for k in 1:2:(np - 1)
            ap_hijo = hijos[k]
            padres_idx = findall(==(ap_hijo), padres)
            pi = isempty(padres_idx) ? rand(rng_dibujo, 1:np) : rand(rng_dibujo, padres_idx)
            for h in (k, k + 1)
                push!(xc, pi); push!(yc, -(g - 1))
                push!(xc, h);  push!(yc, -g)
                push!(xc, NaN); push!(yc, NaN)
            end
        end
    end
    plot!(p_arbol, xc, yc; color = :black, alpha = 0.5, lw = 1.1, label = false)

    for ap in 1:10
        xs = Float64[]; ys = Float64[]
        for g in 0:G
            for i in 1:N
                if hist[g + 1][i] == ap
                    push!(xs, i); push!(ys, -g)
                end
            end
        end
        if !isempty(xs)
            scatter!(p_arbol, xs, ys; color = COLORES_DANE[ap], markersize = 8,
                     markerstrokecolor = :black, markerstrokewidth = 1,
                     label = APELLIDOS_DANE_ASCII[ap])
        end
    end
    p_arbol
end

# ╔═╡ e0000004-0004-0004-0004-000000000022
begin
    nom_ganador = APELLIDOS_DANE_ASCII[arbol.ganador]
    ext_ord = sort(collect(arbol.caidos), by = x -> x[2])
    filas_ext = String[]
    for (pos, (ap_id, gen_ext)) in enumerate(ext_ord)
        push!(filas_ext, "$pos. **$(APELLIDOS_DANE_ASCII[ap_id])** (extinto en Gen $gen_ext)")
    end
    resumen_ext = isempty(filas_ext) ? "*Ningún apellido se extinguió aún.*" : join(filas_ext, "<br>")

    md"""
    ### Resultado del árbol

    | Campo | Valor |
    |:--|:--|
    | Estado final | **$(arbol.fijado ? "Fijación de un apellido" : "Sin fijación en $(arbol.generaciones) generaciones")** |
    | Apellido dominante | **$(nom_ganador)** |
    | Generaciones | **$(arbol.generaciones)** |
    | Modo inicial | **$(modo_dane)** |
    | Modo de probabilidad | **$(modo_prob_dane)** |

    **Cronología de extinción:**<br>
    $(resumen_ext)
    """
end

# ╔═╡ Cell order:
# ╟─e0000004-0004-0004-0004-000000000001
# ╟─e0000004-0004-0004-0004-000000000002
# ╟─e0000004-0004-0004-0004-000000000010
# ╟─e0000004-0004-0004-0004-000000000011
# ╟─e0000004-0004-0004-0004-000000000012
# ╟─e0000004-0004-0004-0004-000000000013
# ╟─e0000004-0004-0004-0004-000000000014
# ╟─e0000004-0004-0004-0004-000000000020
# ╠═e0000004-0004-0004-0004-000000000021
# ╟─e0000004-0004-0004-0004-000000000022
