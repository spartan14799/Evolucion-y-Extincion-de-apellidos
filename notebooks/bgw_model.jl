### A Pluto.jl notebook ###
# v1.0.3

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

	
	# Constantes y utilidades compartidas por todo el cuaderno

	# Duración en años de una generación supuesta
	anios_por_generacion = 28

	# Fuente en negrita para las anotaciones dentro de las graficas.
	fuente_negrita = Sys.iswindows() ? "arialbd" : "sans-serif"

	# Formato de miles con punto
	fmt_miles(n::Integer) = replace(string(n), r"(?<=\d)(?=(\d{3})+$)" => ".")

	# Clasifica el régimen según la media de descendencia m = Σ k·pₖ.
	# La tolerancia es simétrica alrededor de m = 1 para absorber el error de
	# redondeo en punto flotante de la suma que calcula funcion: "esperanza".
	function clasificar_regimen(m::Real; tol::Real = 1e-6)
		if isapprox(m, 1.0; atol = tol)
			return "Crítico"
		elseif m < 1.0
			return "Subcrítico"
		else
			return "Supercrítico"
		end
	end

	# Probabilidad exacta de extinción hasta la generación n, para n = 0, 1, …, n_max:
	#
	#     P(Z_n = 0 | Z_0 = z0) = [ g_n(0) ]^z0 ,   con   g_{n+1}(0) = g( g_n(0) ),  g_0(0) = 0.
	#
	# (sección 4.3 del encabezado). Se itera la FGP del proyecto, `fgp`, partiendo de 0.
	# La potencia z0 se evalúa como exp(z0·ln s) para no acumular error de redondeo
	# cuando z0 es grande (del orden de 10⁵ en el caso de los apellidos).
	function prob_extincion_curva(modelo::ModeloPBGW, n_max::Int, z0::Int)
		s = 0.0                                   # s = g_n(0); empieza en g_0(0) = 0
		curva = Vector{Float64}(undef, n_max + 1)
		for n in 0:n_max
			# min(1, ·) protege frente a un redondeo que dejara s ligeramente por encima de 1
			curva[n + 1] = s <= 0.0 ? 0.0 : min(1.0, exp(z0 * log(s)))
			s = fgp(modelo, s)                    # g_{n+1}(0) = g( g_n(0) )
		end
		return curva
	end

	# P(Z_n = 0 | Z_0 = z0) en una sola generación n.
	prob_extincion_hasta(modelo::ModeloPBGW, n::Int, z0::Int) =
		prob_extincion_curva(modelo, n, z0)[end]

	# Cuantas generaciones tiene que pasar para haya al menos una probabilidad del
	# 50% de que ya se haya extinguido.
	# Devuelve `nothing` si no se alcanza en n_max generaciones (p. ej. si la
	# extinción no es segura, o si m está muy cerca de 1 y z0 es muy grande).
	function mediana_extincion_teorica(modelo::ModeloPBGW, z0::Int; n_max::Int = 20_000)
		s = 0.0
		for n in 1:n_max
			s = fgp(modelo, s)                    # s = g_n(0)
			# [g_n(0)]^z0 ≥ 1/2  ⇔  z0·ln g_n(0) ≥ ln(1/2)
			if s > 0.0 && z0 * log(s) >= log(0.5)
				return n
			end
		end
		return nothing
	end
end

# ╔═╡ 00000002-0002-0002-0002-000000000085
TableOfContents(title = "Índice", indent = true, depth = 4, aside = true)

# ╔═╡ 00000002-0002-0002-0002-000000000002
md"""
# Modelo de Bienaymé–Galton–Watson (PBGW)
## Extinción y supervivencia de apellidos como cadena de Markov de ramificación

---

Este cuaderno estudia el **proceso de ramificación de Bienaymé–Galton–Watson** (PBGW) como una **cadena de Markov en tiempo discreto** sobre los enteros no negativos, y lo aplica al problema que motivó su introducción en el siglo XIX: determinar si un linaje que se transmite por línea paterna, como un apellido, se **extingue** o **persiste** en el tiempo. Bienaymé (1845) y, de forma independiente, Galton y Watson (1875) formularon el modelo para analizar la desaparición de apellidos aristocráticos.

El cuaderno combina tres elementos:

1. **Fundamento teórico.** Definición del proceso, su caracterización como cadena de Markov (espacio de estados, distribución inicial y matriz de transición $P$), la función generadora de probabilidad (FGP) como herramienta de cálculo y el teorema de extinción.
2. **Verificación computacional.** Cada resultado teórico se contrasta con las funciones del proyecto (`src/markovModels.jl` y `src/numericalMethods.jl`) y con simulación Monte Carlo.
3. **Aplicación a los apellidos más frecuentes de Colombia**, con la matriz de transición construida de forma explícita.

### Contenido de la parte teórica

| Sección | Tema |
|:--|:--|
| 1 | Definición del proceso y notación |
| 2 | El PBGW como cadena de Markov: estados, distribución inicial y matriz $P$ (con visualización interactiva) |
| 3 | La función generadora de probabilidad |
| 4 | Iteración de la FGP y ley de $Z_n$ |
| 5 | Probabilidad de extinción y regímenes |
| 6 | Clasificación de estados y comportamiento asintótico |
| 7 | Truncamiento del espacio de estados y métodos numéricos |
| 8 | Supuestos y limitaciones del modelo aplicado a apellidos |
| 9 | Referencias |
"""

# ╔═╡ 00000002-0002-0002-0002-000000000074
md"""
---

## Parte I. Fundamento teórico

### 1. Definición del proceso

Sea $\{X_{n,i}\}_{n \ge 0,\; i \ge 1}$ una familia de variables aleatorias independientes e idénticamente distribuidas (i.i.d.) con valores en $\mathbb{N} = \{0, 1, 2, \dots\}$ y ley común

$$p_k = \mathbb{P}(X = k), \qquad k \ge 0, \qquad \sum_{k \ge 0} p_k = 1 .$$

A $\{p_k\}$ se le llama **ley de descendencia**. La variable $X_{n,i}$ es el número de hijos del individuo $i$ de la generación $n$. El **proceso de Bienaymé–Galton–Watson** es la sucesión $\{Z_n\}_{n \ge 0}$ definida por recurrencia:

$$Z_0 = z_0, \qquad Z_{n+1} = \sum_{i=1}^{Z_n} X_{n,i}, \qquad n \ge 0,$$

con la convención de que una suma vacía vale $0$. $Z_n$ es el número de individuos de la generación $n$; para nuestro caso, el número de portadores del apellido.

**Observación (qué es y qué no es la cadena de Markov).** El objeto a analizar como cadena de Markov no es el árbol genealógico completo, sino la sucesión de **tamaños generacionales** $\{Z_n\}$. El árbol contiene más información (quién es hijo de quién), que es irrelevante para la pregunta de extinción.

**Hipótesis.**

* (H1) Las variables $X_{n,i}$ son independientes entre individuos y entre generaciones: lo que le ocurre a un individuo no afecta la descendencia de los demás.
* (H2) Todas tienen la misma ley $\{p_k\}$: La cadena es homogénea.
* (H3) La media $m = \mathbb{E}[X] = \sum_k k\,p_k$ es finita.
* (H4, no degeneración) $p_1 < 1$. Si $p_1 = 1$, cada individuo tiene exactamente un hijo, $Z_n \equiv z_0$ y no hay extinción.

**Notación y correspondencia con el código.**

| Símbolo | Significado | En el código |
|:--|:--|:--|
| $Z_n$ | tamaño de la generación $n$ | vector devuelto por `simular_generaciones` (posición $n+1$) |
| $X$ | número de hijos de un individuo | — |
| $p_k = \mathbb{P}(X=k)$ | ley de descendencia | `modelo.pk[k+1]` (los vectores de Julia empiezan en 1) |
| $m = \mathbb{E}[X]$ | media de descendencia | `esperanza(modelo)` |
| $\sigma^2 = \mathrm{Var}(X)$ | varianza de la descendencia | — |
| $g(t) = \mathbb{E}[t^X]$ | función generadora de probabilidad | `fgp(modelo, t)` |
| $q$ | probabilidad de extinción partiendo de un individuo | `probabilidad_extincion(modelo)` |
| $T = \inf\{n \ge 0 : Z_n = 0\}$ | tiempo de extinción | primera posición con valor $0$, menos uno |
"""

# ╔═╡ 00000002-0002-0002-0002-000000000075
md"""
### 2. El PBGW como cadena de Markov

Una cadena de Markov homogénea en tiempo discreto queda completamente caracterizada por tres objetos: el **espacio de estados**, la **distribución inicial** y la **matriz de transición**.

#### 2.1 Espacio de estados

$$S = \mathbb{N} = \{0, 1, 2, \dots\},$$

donde el estado $i$ significa "la generación actual tiene exactamente $i$ individuos".

#### 2.2 Distribución inicial

El proceso comienza con $z_0$ individuos de forma determinista, así que la distribución inicial es concentrar toda la probabilidad en $z_0$:

$$\mu_0 = \delta_{z_0}, \qquad \text{es decir,} \qquad \mathbb{P}(Z_0 = z_0) = 1 .$$

Lo habitual es $z_0 = 1$ (un único ancestro). En el cuaderno, $\mu_0$ se representa como un vector fila $\pi_0$ con un $1$ en la posición $z_0$ y ceros en las demás.

#### 2.3 Matriz de transición

Si $Z_n = i$, entonces $Z_{n+1}$ es la suma de $i$ copias independientes de $X$. Por tanto

$$P(i,j) = \mathbb{P}(Z_{n+1} = j \mid Z_n = i) = p_j^{*i} = \sum_{\substack{k_1, \dots, k_i \ge 0 \\ k_1 + \cdots + k_i = j}} p_{k_1} p_{k_2} \cdots p_{k_i},$$

donde $p^{*i}$ es la **convolución $i$-ésima** de la ley $p$ consigo misma, con la convención $p^{*0} = \delta_0$. La convolución es simplemente **contar combinaciones**: $p_j^{*i}$ es la probabilidad de que $i$ individuos independientes sumen exactamente $j$ hijos entre todos.

La estructura por filas es:

* **Fila $0$.** $P(0,0) = 1$ y $P(0,j) = 0$ para $j \ge 1$. Sin individuos no hay quien se reproduzca: el $0$ es un estado **absorbente** (extinción).
* **Fila $1$.** $P(1,j) = p_j$. La primera fila con individuos es exactamente la ley de descendencia.
* **Fila $i \ge 2$.** $P(i,j) = p_j^{*i}$.

$$P = \begin{pmatrix}
1 & 0 & 0 & 0 & \cdots \\
p_0 & p_1 & p_2 & p_3 & \cdots \\
p_0^2 & 2p_0p_1 & p_1^2 + 2p_0p_2 & 2p_0p_3 + 2p_1p_2 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \ddots
\end{pmatrix}$$

Cada fila suma $1$; la fila $i$ tiene media $i\,m$ y varianza $i\,\sigma^2$, por ser la distribución de una suma de $i$ variables i.i.d.

#### 2.4 Por qué es una cadena de Markov

El vector $\left(Z_0, \dots, Z_n\right)$ es función exclusivamente de las variables $\{X_{k,l} : k < n\}$, que son independientes de $\{X_{n,l}\}_{l \ge 1}$. Dado $Z_n = i$, el valor de $Z_{n+1}$ depende solo de $X_{n,1}, \dots, X_{n,i}$. Entonces

$$\mathbb{P}(Z_{n+1} = j \mid Z_n = i, Z_{n-1} = i_{n-1}, \dots, Z_0 = i_0) = \mathbb{P}\Big(\sum_{l=1}^{i} X_{n,l} = j\Big) = p_j^{*i} = P(i,j).$$

Esta es la propiedad de Markov. Como el lado derecho no depende de $n$, la cadena es **homogénea en el tiempo**.

#### 2.5 Distribución en la generación $n$ y propiedad de ramificación

Si $\pi_n(j) = \mathbb{P}(Z_n = j)$ y $\pi_n$ se escribe como vector fila, las ecuaciones de Chapman–Kolmogorov dan

$$\pi_{n+1} = \pi_n P, \qquad \pi_n = \pi_0 P^n .$$

Además, la cadena tiene la **propiedad de ramificación**: partiendo de $i$ individuos, el proceso es la suma de $i$ procesos independientes, cada uno iniciado en un individuo. En particular, $\mathbb{P}_i(Z_n = 0) = \mathbb{P}_1(Z_n = 0)^i$. Esta propiedad justifica que la extinción con $z_0$ individuos iniciales tenga probabilidad $q^{z_0}$ (sección 5).
"""

# ╔═╡ 00000002-0002-0002-0002-000000000084
md"""
#### 2.6 Visualización interactiva de la matriz $P$ y de la cadena

Las celdas siguientes muestran de forma explícita la matriz de transición truncada (sección 7.1), el grafo de estados y la evolución de la distribución $\pi_n = \pi_0 P^n$ (sección 2.5). Se controlan con el tope de estados $M$, el estado inicial $Z_0$ y la generación $n$. La tasa de natalidad $\lambda$ es la del selector de apellidos de la Parte III.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000062
md"**Tope de estados a visualizar ($M$):** $(@bind M_slider Slider(4:1:12, default=8, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000063
md"**Estado inicial $Z_0$ para la cadena:** $(@bind z0_markov_slider Slider(1:1:8, default=2, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000083
begin
	# Matriz de transición truncada del PBGW con descendencia Poisson(λ).
	#
	# Estados 0, 1, …, M (posiciones 1, …, M+1 de la matriz); el último es "M+".
	#   Fila 0:  P(0,0) = 1  (absorbente).
	#   Fila i:  con i portadores la descendencia total es Poisson(i·λ) (propiedad
	#            reproductiva de Poisson), así que P(i,j) = P(Poisson(iλ) = j) para
	#            j = 0, …, M-1, y P(i, M+) = P(Poisson(iλ) ≥ M) es la cola agregada.
	#
	# `max_hijos` debe ser grande para que la cola quede bien representada: el
	# constructor ModeloPBGW corta la ley en max_hijos y RENORMALIZA, de modo que un
	# corte en M distorsionaría tanto P(i,j) como la columna M+ (sección 7.2).
	function matriz_transicion_truncada(lambda::Float64, M::Int; max_hijos::Int = 200)
		P = zeros(Float64, M + 1, M + 1)
		P[1, 1] = 1.0                                   # estado 0 absorbente
		for i in 1:M
			pk_i = ModeloPBGW(i * lambda; max_hijos = max_hijos).pk   # pk_i[k+1] = P(Poisson(iλ) = k)
			P[i + 1, 1:M] = pk_i[1:M]                   # P(i, j) para j = 0, …, M-1
			P[i + 1, M + 1] = sum(pk_i[(M + 1):end])    # P(i, ≥ M): toda la cola va a "M+"
		end
		return P
	end
end

# ╔═╡ 00000002-0002-0002-0002-000000000072
md"""
##### Diagrama de estados de la cadena (grafo de transición)

La cadena se visualiza como una **red dirigida** (ejemplo con $M = 3$: estados $0, 1, 2, 3+$):

- **Nodos.** Estados de la cadena (cantidad de portadores vivos). El estado **0** (rojo) es el pozo absorbente.
- **Flechas dirigidas.** Saltos de transición en un paso (una generación), rotulados con su probabilidad $P_{ij}$. Solo se dibujan las transiciones con probabilidad superior a 0.025.
- **Flechas inferiores rojas y naranjas.** Camino de **decaimiento** ($i \to j$ con $j < i$); la roja es la que conduce a la extinción. Son las más gruesas en régimen subcrítico.
- **Flechas superiores azules.** Camino de **crecimiento** del linaje ($i \to j$ con $j > i$).
- **Autolazos.** Transición $i \to i$ (el tamaño no cambia). En el estado 0 tiene probabilidad $1.00$ y no hay flecha de salida.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000066
md"**Paso de tiempo / generación ($n$) para el vector de estado $\boldsymbol{\pi}_n$:** $(@bind gen_pi_slider Slider(0:1:25, default=3, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000076
md"""
### 3. La función generadora de probabilidad (FGP)

#### 3.1 Definición

Para una variable aleatoria $X$ con valores en $\mathbb{N}$ y $p_k = \mathbb{P}(X = k)$, su **función generadora de probabilidad** es

$$g(t) = \mathbb{E}[t^X] = \sum_{k=0}^{\infty} p_k\, t^k, \qquad t \in [0,1].$$

Es una serie de potencias cuyos **coeficientes son las probabilidades**: contiene exactamente la misma información que la sucesión $\left(p_0, p_1, p_2, \dots\right)$, escrita como un polinomio. La variable $t$ no representa ningún valor aleatorio: solo marca la posición $k$ de cada coeficiente.

#### 3.2 Propiedades que se leen directamente

| Propiedad | Consecuencia |
|:--|:--|
| $g(0) = p_0$ | Probabilidad de que un individuo no tenga hijos |
| $g(1) = 1$ | Las probabilidades suman $1$ (siempre sirve como verificación) |
| $g'(1) = \sum_k k\,p_k = m$ | Media de descendencia |
| $g$ es creciente y **convexa** en $[0,1]$ | $g''(t) = \sum_k k(k-1)p_k t^{k-2} \ge 0$; es estrictamente convexa si $p_0 + p_1 < 1$ |

#### 3.3 La propiedad clave: las sumas independientes se convierten en productos

Si $X$ e $Y$ son independientes,

$$\mathbb{E}\big[t^{X+Y}\big] = \mathbb{E}\big[t^{X}\, t^{Y}\big] = \mathbb{E}\big[t^{X}\big]\,\mathbb{E}\big[t^{Y}\big].$$

La FGP de la suma es el producto de las FGP. Al multiplicar dos polinomios se multiplica cada término del primero por cada término del segundo y los **exponentes se suman**: el exponente $a + b$ representa "$a$ hijos del primer individuo más $b$ hijos del segundo". Al agrupar términos semejantes, el coeficiente de $t^j$ reúne (con sus probabilidades) todas las parejas que suman $j$. Por eso multiplicar polinomios equivale a convolucionar leyes.

#### 3.4 La fila $i$ de $P$ sale de $g(t)^i$

Dado $Z_n = i$, la siguiente generación es la suma de $i$ copias independientes de $X$, y su FGP es $g(t) \cdots g(t) = g(t)^i$. Por tanto

$$\sum_{j \ge 0} P(i,j)\, t^j = \big[g(t)\big]^{i}.$$

**Receta para obtener la fila $i$:** (1) escribir $g(t) = \sum_k p_k t^k$; (2) calcular el polinomio $[g(t)]^i$; (3) leer como $P(i,j)$ el coeficiente de $t^j$. Si $X$ pertenece a una familia cerrada bajo sumas (Bernoulli, Poisson, binomial negativa), el paso (2) puede omitirse: la suma de $i$ copias pertenece a la misma familia.

#### 3.5 Ejemplo resuelto: descendencia uniforme en $\{0, 1, 2\}$

Aquí $p_0 = p_1 = p_2 = 1/3$ y $g(t) = \tfrac13 (1 + t + t^2)$. La media es $m = 1$ (caso crítico). Los polinomios de las primeras filas son:

| $i$ | $\big(1 + t + t^2\big)^i$ | Se divide entre |
|:--|:--|:--|
| 1 | $1 + t + t^2$ | $3$ |
| 2 | $1 + 2t + 3t^2 + 2t^3 + t^4$ | $9$ |
| 3 | $1 + 3t + 6t^2 + 7t^3 + 6t^4 + 3t^5 + t^6$ | $27$ |

Para $i = 2$ hay $3 \times 3 = 9$ combinaciones equiprobables $\left(a,b\right)$ de hijos de cada individuo. Las parejas cuyo total es $0, 1, 2, 3, 4$ son $1, 2, 3, 2, 1$ respectivamente; de ahí los coeficientes de la segunda fila. Si se agrupa en un último estado "$3$ o más" todo lo que supera $2$ (truncamiento con $M = 3$), la matriz de transición es

$$P = \begin{pmatrix}
1 & 0 & 0 & 0 \\
\tfrac{1}{3} & \tfrac{1}{3} & \tfrac{1}{3} & 0 \\
\tfrac{1}{9} & \tfrac{2}{9} & \tfrac{3}{9} & \tfrac{3}{9} \\
\tfrac{1}{27} & \tfrac{3}{27} & \tfrac{6}{27} & \tfrac{17}{27}
\end{pmatrix},$$

con $\mu_0 = (0, 1, 0, 0)$ si se parte de un individuo. La última entrada de la fila 2 es $\tfrac{2+1}{9}$ (totales $3$ y $4$) y la de la fila 3 es $\tfrac{7+6+3+1}{27}$ (totales $3$ a $6$).

#### 3.6 Dos usos distintos de $g$: potencia frente a composición

* **Potencia $g^i$** (multiplicar): más individuos **en la misma generación**. Da la fila $i$ de $P$.
* **Composición $g \circ g$**: más **generaciones**. Da la ley de $Z_2$ partiendo de un individuo (sección 4).
"""

# ╔═╡ 00000002-0002-0002-0002-000000000077
md"""
### 4. Iteración de la FGP y ley de $Z_n$

#### 4.1 Teorema (iteración)

Sea $Z_0 = 1$ y $g_n(s) = \mathbb{E}[s^{Z_n}]$. Entonces

$$g_1 = g, \qquad g_{n+1} = g_n \circ g = g \circ g_n, \qquad \text{es decir,} \qquad g_n = \underbrace{g \circ g \circ \cdots \circ g}_{n \text{ veces}} .$$

*Demostración.* Se pueden condicionar dos instantes distintos.

* Condicionando en $Z_n = i$, la FGP de $Z_{n+1}$ es $g(s)^i$ (sección 3.5). Luego $g_{n+1}(s) = \sum_i \mathbb{P}(Z_n = i)\, g(s)^i = g_n\big(g(s)\big)$, pues $\sum_i \mathbb{P}(Z_n = i)\,t^i$ es, por definición, $g_n(t)$ evaluada en $t = g(s)$.
* Condicionando en $Z_1 = i$ y usando la propiedad de ramificación, $Z_{n+1}$ es la suma de $i$ copias independientes de $Z_n$ iniciadas en un individuo, con FGP $g_n(s)^i$. Luego $g_{n+1}(s) = \sum_i p_i\, g_n(s)^i = g\big(g_n(s)\big)$. $\square$

Si $Z_0 = z_0$, por la propiedad de ramificación $\mathbb{E}[s^{Z_n}] = g_n(s)^{z_0}$.

#### 4.2 Momentos

Derivando $g_n(s) = g\big(g_{n-1}(s)\big)$ en $s = 1$ y usando $g_{n-1}(1) = 1$ resulta $g_n'(1) = g'(1)\, g_{n-1}'(1) = m\, g_{n-1}'(1)$. De aquí,

$$\mathbb{E}[Z_n] = z_0\, m^{n}, \qquad \mathrm{Var}(Z_n) = \begin{cases} z_0\,\sigma^2\, m^{n-1}\,\dfrac{m^{n} - 1}{m - 1}, & m \ne 1, \\[2mm] z_0\,\sigma^2\, n, & m = 1. \end{cases}$$

La media decae o crece de forma geométrica con razón $m$; es lo que representa la curva discontinua $\mathbb{E}[Z_n] = Z_0\,\lambda^n$ en la gráfica de trayectorias en escala logarítmica.

#### 4.3 Distribución exacta del tiempo de extinción

Como $0$ es absorbente, $\{T \le n\} = \{Z_n = 0\}$ y, por la propiedad de ramificación,

$$\mathbb{P}(T \le n) = \mathbb{P}(Z_n = 0) = \big[g_n(0)\big]^{z_0}.$$

La fórmula es **exacta** y se evalúa con $n$ aplicaciones de `fgp` a partir de $0$, mediante $g_{n+1}(0) = g\big(g_n(0)\big)$. El cuaderno la usa para contrastar la simulación Monte Carlo con la teoría.

#### 4.4 Orden de magnitud del tiempo de extinción (caso subcrítico)

Si $m < 1$ y $\mathbb{E}[X \ln X] < \infty$ (condición que cumple la ley de Poisson), existe una constante $c > 0$ tal que $1 - g_n(0) \sim c\, m^{n}$ cuando $n \to \infty$ (Athreya y Ney, 1972). Entonces

$$\mathbb{P}(T \le n) = \big(1 - (1 - g_n(0))\big)^{z_0} \approx \exp\!\big(-z_0\, c\, m^{n}\big),$$

y la mediana de $T$ cumple $z_0\, c\, m^{n} \approx \ln 2$, es decir,

$$T_{\text{med}} \approx \frac{\ln z_0}{\ln(1/m)} + \frac{\ln(c/\ln 2)}{\ln(1/m)} .$$

El término dominante, $\ln z_0 / \ln(1/m)$, es el que el cuaderno usa como "aproximación determinista"; la constante $c$ depende de la ley y aporta un desplazamiento de orden uno. **El tiempo de extinción crece con el logaritmo de la población inicial:** multiplicar $z_0$ por $k$ añade aproximadamente $\ln k / \ln(1/m)$ generaciones.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000003
md"""
### 5. Probabilidad de extinción y regímenes

#### 5.1 Definición

$$q = \mathbb{P}(\text{extinción} \mid Z_0 = 1) = \mathbb{P}(T < \infty) = \lim_{n \to \infty} \mathbb{P}(Z_n = 0) = \lim_{n \to \infty} g_n(0).$$

El límite existe porque los sucesos $\{Z_n = 0\}$ son crecientes: el $0$ es absorbente.

#### 5.2 Teorema (Galton–Watson)

Sea $p_1 < 1$. La probabilidad de extinción $q$ es la **menor raíz no negativa de la ecuación**

$$g(s) = s, \qquad s \in [0, 1].$$

Además, $q = 1$ si y solo si $m \le 1$. Si $m > 1$, entonces $q < 1$ y $q$ es la única raíz de $g(s) = s$ en $[0, 1)$.

#### 5.3 Demostración

**(i) $q$ es un punto fijo.** Se condiciona en el número de hijos del ancestro y se usa la propiedad de ramificación (cada hijo inicia un linaje independiente que se extingue con probabilidad $q$):

$$q = \sum_{k \ge 0} p_k\, q^{k} = g(q).$$

**(ii) $q$ es el menor punto fijo.** Sea $q_n = g_n(0) = \mathbb{P}(Z_n = 0)$. Se cumple $q_0 = 0$ y $q_{n+1} = g(q_n)$. Si $s^\ast \in [0,1]$ es cualquier punto fijo de $g$, como $g$ es creciente, $q_n \le s^\ast$ implica $q_{n+1} = g(q_n) \le g(s^\ast) = s^\ast$. Por inducción, $q_n \le s^\ast$ para todo $n$ y, pasando al límite, $q \le s^\ast$.

**(iii) Localización de las raíces.** Sea $h(s) = g(s) - s$. Es convexa, con $h(0) = p_0 \ge 0$, $h(1) = 0$ y $h'(1) = m - 1$.

* Si $m \le 1$: $h'(1) \le 0$ y, por convexidad (estricta bajo la hipótesis $p_1 < 1$), $h > 0$ en $[0,1)$. La única raíz es $1$, luego $q = 1$.
* Si $m > 1$: $h'(1) > 0$, así que $h < 0$ a la izquierda de $1$, mientras que $h(0) = p_0 \ge 0$. Por el teorema del valor intermedio hay una raíz en $[0,1)$, única por convexidad. $\square$

**Caso $Z_0 = z_0$.** Los $z_0$ linajes son independientes y todos deben extinguirse:

$$\mathbb{P}(\text{extinción} \mid Z_0 = z_0) = q^{z_0}.$$

Por eso el cuaderno compara la frecuencia observada con $q^{Z_0}$, y no con $q$, cuando $Z_0 > 1$.

#### 5.4 Los tres regímenes

| Régimen | Condición | $q$ | $\mathbb{E}[Z_n] = z_0 m^n$ | Función del proyecto |
|:--|:--|:--|:--|:--|
| **Subcrítico** | $m < 1$ | $q = 1$ | decae geométricamente | `es_extincion_segura() → true` |
| **Crítico** | $m = 1$ | $q = 1$ | constante | `es_extincion_segura() → true` |
| **Supercrítico** | $m > 1$ | $q < 1$ | crece geométricamente | `es_extincion_segura() → false` |

* **Subcrítico.** El linaje se contrae y desaparece con rapidez geométrica de razón $m$.
* **Crítico.** La extinción es segura aunque la media se mantenga constante: las trayectorias que sobreviven alcanzan tamaños cada vez mayores que compensan a las extintas. Si $\sigma^2 < \infty$, $n\,\mathbb{P}(Z_n > 0) \to 2/\sigma^2$ (Kolmogorov), de modo que la extinción es muy lenta.
* **Supercrítico.** Con probabilidad $1 - q$ el linaje crece sin cota (dicotomía de la sección 6).
* **Caso degenerado.** Si $p_1 = 1$ (cada individuo tiene un único hijo), $m = 1$ pero $q = 0$: no hay extinción. La hipótesis $p_1 < 1$ del teorema excluye este caso; el cuaderno lo trata de forma explícita.

**Comportamiento cerca del umbral.** Al aproximarse a $m = 1$ por la derecha, $q \approx 1 - 2(m-1)/\sigma^2$. La función $q(m)$ es **continua pero no derivable en $m = 1$**: ese cambio cualitativo es la *transición de fase* que muestra la gráfica de $q$ frente a $m$.

#### 5.5 Cálculo numérico de $q$

Para $m > 1$ se resuelve $h(t) = g(t) - t = 0$. Como $t = 1$ siempre es raíz, se busca en $[0,\, 1 - \varepsilon]$ con $\varepsilon = 10^{-6}$. En ese intervalo $h(0) = p_0 > 0$ y $h(1-\varepsilon) < 0$, de modo que hay **cambio de signo** y el método de bisección converge a $q$. Para $m \le 1$, `probabilidad_extincion` devuelve directamente $1$. Una alternativa equivalente es la iteración de punto fijo $q_{n+1} = g(q_n)$ desde $q_0 = 0$, que converge de forma monótona a $q$ (paso (ii) de la demostración).
"""

# ╔═╡ 00000002-0002-0002-0002-000000000004
md"""
### 6. Clasificación de estados y comportamiento asintótico

Se supone $p_0 > 0$ (con descendencia de Poisson, $p_0 = e^{-\lambda} > 0$ siempre).

| Propiedad | Resultado | Justificación |
|:--|:--|:--|
| Estado $0$ | **absorbente** y recurrente | $P(0,0) = 1$ |
| Estados $i \ge 1$ | **transitorios** | desde $i$ se llega a $0$ con probabilidad $p_0^{\,i} > 0$ y de $0$ no se regresa |
| Irreducibilidad | **no irreducible** | de los estados positivos se llega a $0$, pero no al revés |
| Clases de comunicación (Poisson) | $\{0\}$ cerrada y $\{1, 2, 3, \dots\}$ abierta | como $p_k > 0$ para todo $k$, los estados positivos se comunican entre sí |
| Periodicidad | **aperiódica** | $P(i,i) > 0$ para todo $i$ |
| Distribución estacionaria | única: $\delta_0$ | toda masa en estados transitorios termina disipándose |

**Dicotomía fundamental.** Si $p_1 < 1$ y $Z_0 = 1$, entonces $\mathbb{P}(Z_n \to 0) = q$ y $\mathbb{P}(Z_n \to \infty) = 1 - q$. Ningún estado positivo es visitado infinitas veces, así que el proceso no puede estabilizarse en un tamaño finito positivo.

**Límite de la distribución $\pi_n$.**

* Si $m \le 1$: $\pi_n \to \delta_0$ para cualquier distribución inicial. Toda la masa de probabilidad termina concentrada en el estado $0$, y el tiempo de absorción es finito con probabilidad $1$.
* Si $m > 1$: $\pi_n(0) \to q^{z_0}$ y $\pi_n(j) \to 0$ para todo $j \ge 1$. La masa restante, $1 - q^{z_0}$, escapa a infinito: no existe distribución límite sobre $S$.

**Consecuencia para el muestreo.** Como el único estado recurrente es $0$ y no hay distribución estacionaria no trivial, no corresponde usar métodos MCMC (cuyo fin es muestrear una distribución estacionaria o *a posteriori*). Lo que interesa es la **trayectoria dinámica** hacia la extinción, y se estudia con **Monte Carlo directo** (hacia adelante) sobre la cadena $\{Z_n\}$.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000078
md"""
### 7. Truncamiento del espacio de estados y métodos numéricos

#### 7.1 Matriz de transición truncada

La matriz $P$ es infinita. Para mostrarla de forma explícita se trunca a $S_M = \{0, 1, \dots, M\}$, agrupando en el último estado a todos los linajes con $M$ o más portadores:

$$P_M(i,j) = P(i,j) \quad (0 \le j < M), \qquad P_M(i,M) = \sum_{j \ge M} P(i,j) = 1 - \sum_{j < M} P(i,j).$$

Para el estado agregado "$M+$" se usa la fila del estado $M$. Esto es una **aproximación**: el proceso agregado no es exactamente markoviano, porque los estados $M, M+1, \dots$ tienen filas distintas. Tratar "$M+$" como si fuera exactamente $M$ subestima el tamaño del linaje y, por tanto, **sobrestima ligeramente la probabilidad de extinción**. El error es pequeño cuando es poco probable superar $M$ (régimen subcrítico con $M$ moderado); el cuaderno lo cuantifica comparando la curva $\mathbb{P}(Z_n = 0)$ de la cadena truncada con la fórmula exacta $[g_n(0)]^{z_0}$.

#### 7.2 Atajo de Poisson para las filas de $P$

Si $X \sim \mathrm{Poisson}(\lambda)$, la suma de $i$ copias independientes es $\mathrm{Poisson}(i\lambda)$. La fila $i$ de $P$ es entonces

$$P(i,j) = e^{-i\lambda}\,\frac{(i\lambda)^{j}}{j!}, \qquad j \ge 0,$$

sin necesidad de calcular convoluciones. La FGP lo confirma: $[g(t)]^i = \big(e^{\lambda(t-1)}\big)^i = e^{i\lambda(t-1)}$, que es la FGP de una Poisson$\left(i\lambda\right)$.

**Nota de implementación.** El constructor `ModeloPBGW(λ; max_hijos)` corta la ley de Poisson en `max_hijos` y **renormaliza** el vector resultante. Para $\lambda \le 4$ y `max_hijos = 20` el error es inferior a $10^{-8}$. Pero para construir las filas de $P$ hace falta la **cola** $\mathbb{P}(X \ge M)$ de una Poisson$\left(i\lambda\right)$, y para ello hay que cortar muy por encima de $M$ (el cuaderno usa `max_hijos = 200`). Un corte en $M$ seguido de renormalización distorsionaría tanto las probabilidades individuales como la columna "$M+$": por ejemplo, con $\lambda = 0{,}5$ y $M = 3$ la entrada $P(3,3+)$ saldría $0{,}134$ en lugar del valor correcto $0{,}191$.

#### 7.3 Métodos numéricos y simulación

* **Bisección** (`biseccion` en `src/numericalMethods.jl`) sobre $h(t) = g(t) - t$, con tolerancia $10^{-5}$ tanto en el intervalo como en el valor de $h$. Converge linealmente y garantiza una raíz porque hay cambio de signo.
* **Simulación hacia adelante** (`simular_generaciones`): en cada generación se sortea, para cada individuo vivo, su número de hijos con `Categorical(modelo.pk)` y se suman. Es Monte Carlo directo sobre la cadena $\{Z_n\}$.
* **Estimador de extinción** (`frecuencia_extincion`): fracción de $N$ simulaciones con $Z_G = 0$. Estima $\mathbb{P}(Z_G = 0)$ y no $q^{z_0}$: como $\mathbb{P}(Z_G = 0)$ crece con $G$ hacia $q^{z_0}$, el estimador es una cota inferior de la extinción asintótica. Su error estándar es $\sqrt{\hat p(1-\hat p)/N}$.
* **Tope de población** (`tope_poblacion`): en regímenes supercríticos, un linaje que supera el tope se declara establecido (ya no se extinguirá) y se detiene su simulación. Esto acota el costo computacional sin alterar de forma apreciable el resultado de extinción.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000079
md"""
### 8. Supuestos y limitaciones del modelo aplicado a apellidos

La aplicación a apellidos reales exige interpretar con cuidado los parámetros. Deben tenerse presentes las siguientes limitaciones al leer los resultados:

1. **Transmisión patrilineal.** El primer apellido se hereda del padre, de modo que el proceso describe la **línea masculina**: $X$ es el número de **hijos varones** de un hombre y $\lambda$ su media. Si $\mathrm{TGF}$ es la tasa global de fecundidad y la proporción de varones al nacer es cercana a $0{,}512$, entonces $\lambda \approx 0{,}512 \cdot \mathrm{TGF}$, bajo el supuesto de que los hombres tienen en promedio tantos hijos como las mujeres.
2. **Población inicial.** Los conteos usados como $Z_0$ corresponden a **personas** con el apellido, de ambos sexos. Estrictamente, $Z_0$ debería ser el número de **varones** portadores, aproximadamente la mitad. Usar el total sobrestima el tiempo de extinción en $\ln 2 / \ln(1/\lambda)$ generaciones (una generación cuando $\lambda = 0{,}5$).
3. **Generaciones discretas y no solapadas.** El modelo agrupa el tiempo en generaciones sucesivas. La conversión a años usa una duración generacional constante de 28 años, que es un supuesto del modelo y no un dato.
4. **Homogeneidad e independencia.** Todos los portadores comparten la misma ley de descendencia y se reproducen de forma independiente. Se ignoran las diferencias por región, nivel socioeconómico o cohorte, y la tasa se mantiene constante en el tiempo.
5. **Forma de la ley.** La ley de Poisson impone varianza igual a la media. Otras leyes con la misma media dan el mismo $q$ cuando $m \le 1$ (extinción segura), pero distinto tiempo de extinción, porque la constante $c$ de la sección 4.4 depende de toda la ley.
6. **Población cerrada.** No se consideran inmigración, cambios de apellido ni la aparición de nuevos portadores por mecanismos distintos de la descendencia paterna.
7. **Truncamiento numérico.** La ley de Poisson se corta en `max_hijos` y se renormaliza (sección 7.2), con un error despreciable para los valores de $\lambda$ usados.

Estas limitaciones afectan la interpretación cuantitativa. Las conclusiones cualitativas (extinción segura si $m \le 1$ y dependencia logarítmica de $Z_0$) son menos sensibles a ellas.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000080
md"""
### 9. Referencias

- Athreya, K. B. y Ney, P. E. (1972). *Branching Processes*. Springer.
- Bienaymé, I.-J. (1845). De la loi de multiplication et de la durée des familles. *Société Philomathique de Paris, Extraits*, 5, 37–39.
- DANE. *Estadísticas Vitales (EEVV): nacimientos y defunciones*, años 2024 y 2025 (preliminar).
- Galton, F. y Watson, H. W. (1875). On the probability of the extinction of families. *Journal of the Anthropological Institute of Great Britain and Ireland*, 4, 138–144.
- Grimmett, G. R. y Stirzaker, D. R. (2001). *Probability and Random Processes* (3.ª ed.). Oxford University Press. Sección 5.4.
- Harris, T. E. (1963). *The Theory of Branching Processes*. Springer.
- Norris, J. R. (1997). *Markov Chains*. Cambridge University Press.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000005
md"""
---

## Parte II. Análisis interactivo del modelo

### Panel de control

Elija la ley de descendencia y sus parámetros. Todos los cálculos y gráficas siguientes se actualizan automáticamente.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000006
md"**Tipo de distribución de natalidad:** $(@bind dist_tipo Select([\"Poisson\" => \"Poisson (λ)\", \"Binomial\" => \"Binomial (n, p)\"], default=\"Poisson\"))"

# ╔═╡ 00000002-0002-0002-0002-000000000007
md"""
!!! info "Leyes de descendencia disponibles"
    - **Poisson(λ).** Cada portador tiene en promedio λ hijos y no existe un máximo fijo (el cuaderno corta la ley en 20 hijos y renormaliza). Es el modelo clásico y se construye con `ModeloPBGW(λ)`. La media es m = λ y la varianza σ² = λ.

    - **Binomial(n, p).** Cada portador tiene como máximo `n` hijos, cada uno con probabilidad `p`. La media es m = n·p y la varianza σ² = n·p·(1 − p). Resulta apropiada cuando existe un número máximo natural de hijos. Se construye pasando el vector `pk` a `ModeloPBGW(pk)`.
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
	# ── 1. Construcción del modelo de descendencia ──
	if dist_tipo == "Poisson"
		# Poisson(λ) cortada en 20 hijos y renormalizada (sección 7.2 del encabezado)
		modelo = ModeloPBGW(lambda_slider; max_hijos=20)
		param_str = "Poisson(λ=$(lambda_slider))"
	else
		# Binomial(n, p): tiene soporte finito {0, …, n}, así que no requiere truncamiento.
		# pk_binom[k+1] = P(X = k)
		dist_binom = Distributions.Binomial(n_binom_slider, p_binom_slider)
		pk_binom = [pdf(dist_binom, k) for k in 0:n_binom_slider]
		pk_binom = pk_binom ./ sum(pk_binom)  # corrige desvíos de redondeo en la suma
		modelo = ModeloPBGW(pk_binom)
		param_str = "Binomial(n=$(n_binom_slider), p=$(p_binom_slider))"
	end

	# ── 2. Cantidades teóricas (sección 5 del encabezado) ──
	# Media de descendencia m = Σ k·pₖ  (= g'(1))
	m_esperanza = esperanza(modelo)

	# Caso degenerado p₁ = 1 (cada portador tiene exactamente un hijo): entonces
	# Zₙ ≡ Z₀ y nunca hay extinción, aunque m = 1. El teorema de extinción exige
	# p₁ < 1, así que este caso se trata aparte para no reportar q = 1 por error.
	degenerado = length(modelo.pk) >= 2 && isapprox(modelo.pk[2], 1.0; atol = 1e-12)

	# Probabilidad de extinción partiendo de un individuo: menor raíz de g(q) = q
	q_extincion = degenerado ? 0.0 : probabilidad_extincion(modelo)
	extincion_segura = !degenerado && es_extincion_segura(modelo)
	regimen = degenerado ? "Degenerado (p₁ = 1)" : clasificar_regimen(m_esperanza)

	# Con Z₀ portadores iniciales los linajes son independientes y todos deben
	# extinguirse, de modo que P(extinción | Z₀) = q^Z₀  (sección 5.3)
	q_z0 = q_extincion^z0_slider

	texto_segura = if degenerado
		"No (caso degenerado, p₁ = 1)"
	else
		extincion_segura ? "Sí (m ≤ 1)" : "No (m > 1)"
	end

	md"""
	### Análisis del modelo: $(param_str)

	| Parámetro | Valor |
	|:--|:--|
	| Distribución | **$(param_str)** |
	| Media de descendencia (m) | **$(round(m_esperanza, digits=4))** — `esperanza(modelo)` |
	| Régimen | **$(regimen)** |
	| ¿Extinción segura? | **$(texto_segura)** — `es_extincion_segura(modelo)` |
	| Probabilidad de extinción (q) | **$(round(q_extincion, digits=6))** — `probabilidad_extincion(modelo)` |
	| Probabilidad de supervivencia (1 − q) | **$(round(1 - q_extincion, digits=6))** |
	| Población inicial (Z₀) | **$(z0_slider)** |
	| Probabilidad de extinción con Z₀ portadores (q^Z₀) | **$(round(q_z0, digits=6))** |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000015
md"""
### Distribución de descendencia $\{p_k\}$

El vector `modelo.pk` almacena las probabilidades: `pk[1]` es $P(X=0)$, `pk[2]` es $P(X=1)$, y así sucesivamente (los vectores de Julia empiezan en el índice 1). La línea discontinua marca la media $m$.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000016
begin
	# Se recorta la cola de probabilidad despreciable (< 1e-6) para que el eje sea legible
	pk_vec = modelo.pk
	ultimo_significativo = findlast(x -> x > 1e-6, pk_vec)
	k_rango = 0:(ultimo_significativo - 1)        # k = número de hijos
	pk_plot = pk_vec[1:ultimo_significativo]      # pₖ = P(X = k)

	p_dist = bar(
		collect(k_rango), pk_plot,
		xlabel = "Número de hijos (k)",
		ylabel = "Probabilidad P(X = k)",
		title = "Distribución de descendencia — $(param_str)",
		label = false, color = :teal, alpha = 0.8,
		size = (800, 380), xticks = k_rango,
		ylims = (0, 1.2 * maximum(pk_plot))       # margen superior para que ninguna barra toque el borde
	)
	vline!(p_dist, [m_esperanza],
		label = "Media m = $(round(m_esperanza, digits=3))",
		color = :red, linewidth = 2, linestyle = :dash)
	p_dist
end

# ╔═╡ 00000002-0002-0002-0002-000000000017
md"""
### Función generadora de probabilidad $g(t)$ y punto fijo

`probabilidad_extincion()` calcula $q$ resolviendo $g(t) = t$ (sección 5):

- **Panel izquierdo.** $g(t)$ frente a la identidad $y = t$. La probabilidad de extinción es la menor intersección de ambas curvas en $[0,1]$. En régimen subcrítico o crítico la única intersección es $t = 1$ y, por tanto, $q = 1$.
- **Panel derecho.** $h(t) = g(t) - t$. La raíz de $h$ es $q$, hallada con `biseccion()`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000018
begin
	ts = range(0, 1, length=500)
	# g(t) evaluada con la FGP del proyecto, y h(t) = g(t) − t, cuya raíz es q
	g_vals = [fgp(modelo, t) for t in ts]
	h_vals = g_vals .- ts

	# ¿Existe una raíz de g(t) = t distinta de t = 1? (solo si el modelo es supercrítico)
	hay_q_interior = q_extincion < 1.0 - 1e-6
	etiqueta_t1 = hay_q_interior ? "t = 1 (raíz trivial)" : "t = 1 (única raíz, q = 1)"

	p_fgp = plot(layout = (1, 2), size = (950, 440))

	# Panel izquierdo: g(t) frente a la identidad
	plot!(p_fgp[1], ts, g_vals, label = "g(t) = fgp(modelo, t)", color = :blue, linewidth = 2.5)
	plot!(p_fgp[1], ts, ts, label = "y = t", color = :gray, linewidth = 1.5, linestyle = :dash)
	if hay_q_interior
		scatter!(p_fgp[1], [q_extincion], [q_extincion],
			label = "q = $(round(q_extincion, digits=4))",
			color = :red, markersize = 8, markershape = :star5)
	end
	scatter!(p_fgp[1], [1.0], [1.0], label = etiqueta_t1, color = :orange, markersize = 6)
	plot!(p_fgp[1], xlabel = "t", ylabel = "g(t)", title = "FGP y punto fijo",
		legend = :outerbottom)

	# Panel derecho: h(t) = g(t) − t
	plot!(p_fgp[2], ts, h_vals, label = "h(t) = g(t) − t", color = :purple, linewidth = 2.5)
	hline!(p_fgp[2], [0], label = false, color = :gray, linestyle = :dash)
	if hay_q_interior
		scatter!(p_fgp[2], [q_extincion], [0],
			label = "Raíz q = $(round(q_extincion, digits=4))",
			color = :red, markersize = 8, markershape = :star5)
	end
	scatter!(p_fgp[2], [1.0], [0.0], label = etiqueta_t1, color = :orange, markersize = 6)
	plot!(p_fgp[2], xlabel = "t", ylabel = "h(t)", title = "Raíz de h(t) — biseccion()",
		legend = :outerbottom)

	p_fgp
end

# ╔═╡ 00000002-0002-0002-0002-000000000019
md"""
### Trayectorias del linaje — `simular_generaciones()`

Cada trayectoria se genera con `simular_generaciones(modelo, gen; poblacion_inicial=Z₀)`, que simula la cadena $\{Z_n\}$ hacia adelante. Las trayectorias en rojo se extinguieron dentro del horizonte simulado y las verdes siguen vivas al final.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000020
md"**Número de trayectorias:** $(@bind n_traj_slider Slider(1:1:20, default=8, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000021
begin
	rng_traj = MersenneTwister(seed_slider)   # generador reproducible a partir de la semilla
	tope_sim = 10_000                         # tope de población de la simulación (ver nota abajo)

	# Se simulan todas las trayectorias con el mismo generador y en orden
	trayectorias = [simular_generaciones(modelo, gen_slider;
			poblacion_inicial = z0_slider, tope_poblacion = tope_sim, rng = rng_traj)
		for _ in 1:n_traj_slider]

	# Una trayectoria está extinta si termina en 0 al final del horizonte
	extinta_flag = [tr[end] == 0 for tr in trayectorias]
	n_extintas = count(extinta_flag)
	n_vivas = n_traj_slider - n_extintas

	# El eje vertical reserva espacio arriba para la leyenda, de modo que no tape las curvas
	y_max_traj = maximum(maximum.(trayectorias))

	p_traj = plot(
		xlabel = "Generación (n)", ylabel = "Portadores (Zₙ)",
		title = "Trayectorias de Zₙ — $(param_str), Z₀ = $(z0_slider)",
		size = (850, 450), legend = :topleft,
		ylims = (-0.03 * y_max_traj, 1.4 * y_max_traj),
		grid = true, gridstyle = :dash, gridalpha = 0.3
	)

	# Solo se etiqueta la primera trayectoria de cada tipo para mantener breve la leyenda
	i_primera_extinta = findfirst(extinta_flag)
	i_primera_viva = findfirst(!, extinta_flag)
	for (i, tr) in enumerate(trayectorias)
		extinta = extinta_flag[i]
		etiqueta = false
		if extinta && i == i_primera_extinta
			etiqueta = "Extinguida (Zₙ = 0)"
		elseif !extinta && i == i_primera_viva
			etiqueta = "Sobrevive"
		end
		plot!(p_traj, 0:(length(tr) - 1), tr, label = etiqueta,
			color = extinta ? :red : :green,
			alpha = extinta ? 0.55 : 0.8, linewidth = extinta ? 1.5 : 2.0)
	end

	# Nota: cuando un linaje alcanza el tope se declara "establecido" y su población
	# se congela en el tope. La meseta resultante es un artefacto de la simulación,
	# no un comportamiento del modelo, y se señala explícitamente.
	if any(maximum(tr) >= tope_sim for tr in trayectorias)
		hline!(p_traj, [tope_sim],
			label = "Tope de simulación ($(tope_sim)): linaje establecido",
			color = :gray, linestyle = :dot, linewidth = 1.5)
	end

	p_traj
end

# ╔═╡ 00000002-0002-0002-0002-000000000022
begin
	# Probabilidad exacta de que el linaje ya se haya extinguido en la generación
	# gen_slider: [g_G(0)]^Z₀  (sección 4.3). Es la referencia correcta para la
	# frecuencia observada en un horizonte finito.
	p_hasta_gen = prob_extincion_hasta(modelo, gen_slider, z0_slider)
	frec_trayectorias = n_extintas / n_traj_slider

	md"""
	**De $(n_traj_slider) trayectorias:** $(n_extintas) extintas y $(n_vivas) con linaje vivo al final del horizonte.

	| Cantidad | Valor |
	|:--|:--|
	| Frecuencia observada de extinción | **$(round(frec_trayectorias, digits=4))** |
	| Probabilidad exacta de extinción hasta la generación $(gen_slider) | **$(round(p_hasta_gen, digits=4))** |
	| Probabilidad asintótica de extinción, q^Z₀ | **$(round(q_z0, digits=4))** |

	Con pocas trayectorias la frecuencia observada fluctúa mucho; la validación cuantitativa se hace en la sección siguiente con miles de simulaciones.
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000023
md"""
### Validación Monte Carlo — `frecuencia_extincion()`

Se contrasta la estimación empírica (`frecuencia_extincion`, Monte Carlo) con la teoría. Hay que comparar cantidades del mismo tipo:

- `frecuencia_extincion(modelo, G)` estima $\mathbb{P}(Z_G = 0)$, la probabilidad de extinción **hasta la generación $G$**. Su valor exacto es $[g_G(0)]^{Z_0}$ (sección 4.3) y es la referencia contra la que se mide el error de Monte Carlo.
- La probabilidad asintótica $q^{Z_0}$ es el límite de $\mathbb{P}(Z_G = 0)$ cuando $G \to \infty$. La frecuencia empírica nunca la supera y solo se le aproxima si $G$ es grande; cerca del umbral $m = 1$ la convergencia es lenta.

El criterio de concordancia es estadístico: la diferencia entre la frecuencia observada y el valor exacto debe quedar dentro de tres errores estándar, $3\sqrt{p(1-p)/N}$.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000024
md"**Simulaciones Monte Carlo:** $(@bind n_mc_slider Slider(500:500:10000, default=3000, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000025
begin
	rng_mc = MersenneTwister(seed_slider + 7777)

	# Estimación Monte Carlo con la función del proyecto: fracción de simulaciones con Z_G = 0
	freq_mc = frecuencia_extincion(modelo, gen_slider;
		n_sim=n_mc_slider, poblacion_inicial=z0_slider,
		tope_poblacion=10_000, rng=rng_mc)

	# Referencia teórica exacta al mismo horizonte G = gen_slider
	p_exacta = prob_extincion_hasta(modelo, gen_slider, z0_slider)

	# Error estándar del estimador (proporción binomial). El piso 1e-12 evita dividir
	# por cero cuando p_exacta es 0 o 1.
	se_mc = sqrt(max(p_exacta * (1 - p_exacta), 1e-12) / n_mc_slider)
	diferencia = abs(freq_mc - p_exacta)
	z_score = (freq_mc - p_exacta) / se_mc       # desviación en unidades de error estándar

	concordancia = abs(z_score) <= 3 ? "Consistente (a menos de 3 errores estándar)" : "Discrepante (a más de 3 errores estándar)"
	txt_diferencia = "$(round(diferencia, digits=6)) ($(round(z_score, digits=2)) errores estándar)"

	md"""
	| Medida | Valor | Cómo se obtiene |
	|:--|:--|:--|
	| Extinción asintótica, q^Z₀ | **$(round(q_z0, digits=6))** | `probabilidad_extincion(modelo)` elevada a Z₀ |
	| Extinción exacta hasta la generación $(gen_slider) | **$(round(p_exacta, digits=6))** | iteración de `fgp`: [gₙ(0)]^Z₀ |
	| Frecuencia empírica (N = $(n_mc_slider)) | **$(round(freq_mc, digits=6))** | `frecuencia_extincion(modelo, gen_slider)` |
	| Error estándar | **$(round(se_mc, digits=6))** | √(p(1−p)/N) con p exacta |
	| Diferencia respecto a la exacta | **$(txt_diferencia)** | |
	| Concordancia | **$(concordancia)** | |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000026
md"""
### Transición de fase: extinción frente a supervivencia

¿Cómo cambia $q$ al variar la media $m$? En $m = 1$ se produce una **transición de fase**: $q(m)$ es continua pero no derivable (sección 5.4).

Cada punto se calcula con `probabilidad_extincion(ModeloPBGW(λ))`, con descendencia Poisson($\lambda$), para la cual $m = \lambda$:
"""

# ╔═╡ 00000002-0002-0002-0002-000000000027
begin
	lambdas = range(0.1, 4.0, length=100)
	qs = Float64[]
	for lam in lambdas
		# Para cada λ se construye el modelo Poisson y se calcula q con el proyecto:
		# q = 1 si λ ≤ 1; q = menor raíz de g(t) = t (bisección) si λ > 1
		mod_temp = ModeloPBGW(lam; max_hijos=20)
		push!(qs, probabilidad_extincion(mod_temp))
	end

	p_phase = plot(
		lambdas, qs,
		xlabel = "Tasa de natalidad (λ = media m)",
		ylabel = "Probabilidad de extinción (q)",
		title = "Transición de fase: q frente a m",
		label = "q(m) — probabilidad_extincion()",
		color = :darkblue, linewidth = 3,
		size = (800, 420), ylims = (-0.05, 1.1), grid = true
	)
	vline!(p_phase, [1.0], label = "m = 1 (umbral)", color = :red, linestyle = :dash, linewidth = 2)

	# Se marca el modelo actual solo si es Poisson (la curva q(m) es la de Poisson)
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
### Tabla resumen para distintas tasas de natalidad
"""

# ╔═╡ 00000002-0002-0002-0002-000000000029
begin
	lineas_tabla = [
		"| λ (media m) | Régimen | q (extinción) | 1 − q (supervivencia) |",
		"|:--|:--|:--|:--|"
	]
	for lam in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
		mod_t = ModeloPBGW(lam; max_hijos=20)
		m_t = esperanza(mod_t)                  # media de descendencia del modelo
		q_t = probabilidad_extincion(mod_t)     # probabilidad de extinción
		push!(lineas_tabla, "| $(lam) | $(clasificar_regimen(m_t)) | $(round(q_t, digits=4)) | $(round(1 - q_t, digits=4)) |")
	end
	Markdown.parse(join(lineas_tabla, "\n"))
end

# ╔═╡ 00000002-0002-0002-0002-000000000030
md"""
### Probabilidad de absorción con $Z_0 > 1$

Cuando el linaje empieza con $Z_0 > 1$ portadores, cada uno genera un sublinaje **independiente** y todos deben extinguirse para que el apellido desaparezca (propiedad de ramificación, sección 2.5):

$$P(\text{extinción} \mid Z_0) = q^{Z_0}.$$

Esto puede verificarse con `frecuencia_extincion(modelo, gen; poblacion_inicial=Z₀)`.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000031
begin
	z0_vals = 1:10
	q_z0_vals = [q_extincion^z for z in z0_vals]       # q^Z₀ para Z₀ = 1, …, 10

	# Formato de porcentaje que no oculta valores pequeños pero no nulos
	function fmt_pct(x)
		x == 0 && return "0%"
		x < 0.001 && return "<0.1%"
		return string(round(x * 100, digits = 1), "%")
	end

	p_absorcion = bar(
		collect(z0_vals), q_z0_vals,
		xlabel = "Población inicial Z₀",
		ylabel = "P(extinción) = q^Z₀",
		title = "Efecto de Z₀ — $(param_str) (q = $(round(q_extincion, digits=3)))",
		label = false,
		color = [qz < 0.5 ? :green : :red for qz in q_z0_vals],   # verde si la extinción es menos probable que no
		alpha = 0.8, size = (800, 420), xticks = z0_vals,
		ylims = (0, 1.3), yticks = 0:0.25:1                       # margen para etiquetas y leyenda; marcas solo en [0, 1]
	)
	hline!(p_absorcion, [0.5], label = "Umbral 50%", color = :gray, linestyle = :dash)
	for (z, q_z) in zip(z0_vals, q_z0_vals)
		annotate!(p_absorcion, z, q_z + 0.05, text(fmt_pct(q_z), 7, :center))
	end
	p_absorcion
end

# ╔═╡ 00000002-0002-0002-0002-000000000032
md"""
### Resumen y conclusiones de la Parte II

**Funciones del proyecto utilizadas en este cuaderno:**

| Función | Archivo | Propósito |
|:--|:--|:--|
| `ModeloPBGW(λ)` / `ModeloPBGW(pk)` | `src/markovModels.jl` | Construir el modelo con descendencia Poisson o con un vector explícito |
| `esperanza(modelo)` | `src/markovModels.jl` | Media de descendencia $m$ |
| `es_extincion_segura(modelo)` | `src/markovModels.jl` | Verificar si $m \leq 1$ |
| `fgp(modelo, t)` | `src/markovModels.jl` | Función generadora $g(t)$ |
| `probabilidad_extincion(modelo)` | `src/markovModels.jl` | Resolver $g(t) = t$ con bisección |
| `simular_generaciones(modelo, gen)` | `src/markovModels.jl` | Simular una trayectoria de $\{Z_n\}$ |
| `frecuencia_extincion(modelo, gen)` | `src/markovModels.jl` | Monte Carlo: fracción de linajes extintos |
| `biseccion(h, a, b)` | `src/numericalMethods.jl` | Búsqueda de raíces (usada internamente) |
| `secante(h, a, b)` | `src/numericalMethods.jl` | Búsqueda de raíces (alternativa) |

**Conclusiones principales:**

- Todo depende de $m$: si $m \leq 1$ la extinción es segura; si $m > 1$ es posible la supervivencia.
- La FGP y la bisección dan la probabilidad de extinción de forma exacta; la iteración $g_n(0)$ da además la probabilidad de extinción en cualquier generación finita.
- Monte Carlo confirma la teoría: `frecuencia_extincion` coincide con $[g_G(0)]^{Z_0}$ dentro del error estadístico y se aproxima por debajo a `probabilidad_extincion` elevada a $Z_0$.
- Con $Z_0 > 1$ portadores iniciales la probabilidad de extinción disminuye exponencialmente: $q^{Z_0}$.

---
*Cuaderno interactivo del proyecto [Evolución y Extinción de Apellidos](https://github.com/spartan14799/Evolucion-y-Extincion-de-apellidos)*
"""

# ╔═╡ 00000002-0002-0002-0002-000000000040
md"""
---

## Parte III. Aplicación: ¿cuándo se extinguirán los apellidos de Colombia?
### Proceso de Bienaymé–Galton–Watson con los apellidos más frecuentes del país

#### Contexto demográfico

La tasa global de fecundidad (TGF) de Colombia ha descendido de forma sostenida. Según el DANE (Estadísticas Vitales), fue de **1,1 hijos por mujer en 2024** y de **1,0 en 2025** (cifra preliminar), muy por debajo del nivel de reemplazo, de aproximadamente 2,1 hijos por mujer.

Como tradicionalmente el primer apellido lo transmite el padre (supuesto 1 de la sección 8):

- Con una TGF cercana a 1,0 y una proporción de varones al nacer de aproximadamente 0,512, cada hombre que porta el apellido tiene en promedio $\lambda \approx 0{,}51$ hijos varones que continuarán su linaje. El valor por defecto del control, $\lambda = 0{,}5$, resume este caso.
- Como $m = \lambda < 1$, el linaje se encuentra en **régimen subcrítico** y, por el teorema de la sección 5, la extinción matemática es **segura** ($q = 1$).

La pregunta relevante deja de ser *si* el apellido se extingue y pasa a ser **cuándo**: ¿en cuántas generaciones ocurrirá la extinción para los apellidos más numerosos de Colombia? El control de $\lambda$ admite valores mayores que 1 para explorar también el régimen supercrítico.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000041
md"**Seleccione un apellido de Colombia:** $(@bind ap_col_select Select([\"Rodríguez\" => \"1. Rodríguez (754.991 personas)\", \"Martínez\" => \"2. Martínez (602.593 personas)\", \"García\" => \"3. García (574.856 personas)\", \"Gómez\" => \"4. Gómez (569.534 personas)\", \"López\" => \"5. López (565.394 personas)\", \"Hernández\" => \"6. Hernández (502.626 personas)\", \"Pérez\" => \"7. Pérez (493.764 personas)\", \"Sánchez\" => \"8. Sánchez (487.112 personas)\", \"González\" => \"9. González (470.505 personas)\", \"Ramírez\" => \"10. Ramírez (441.179 personas)\"], default=\"Rodríguez\"))"

# ╔═╡ 00000002-0002-0002-0002-000000000042
md"**Tasa de natalidad masculina $\lambda$** (hijos varones promedio por hombre): $(@bind lambda_col_slider Slider(0.1:0.05:1.2, default=0.5, show_value=true))"

# ╔═╡ 00000002-0002-0002-0002-000000000064
begin
	M_dim = M_slider
	lambda_p = lambda_col_slider

	# Matriz de transición explícita (filas = estado actual i, columnas = estado siguiente j)
	P_mat = matriz_transicion_truncada(lambda_p, M_dim)

	# Etiquetas de los estados: "0", "1", …, "M+"
	etiquetas_estados = [string(i) for i in 0:M_dim]
	etiquetas_estados[end] = "$(M_dim)+"

	# Paleta cálida y clara para garantizar máximo contraste con el texto negro
	cmap_mat_markov = cgrad([:white, "#fef08a", "#fdba74", "#f87171", "#ef4444"], [0.0, 0.05, 0.35, 0.7, 1.0])
	fsize_p = M_dim <= 8 ? 10 : 8

	# Los ejes son numéricos (0, …, M) con marcas de texto. Así las coordenadas de las
	# anotaciones coinciden exactamente con el centro de cada celda del mapa de calor.
	p_heatmap_P = heatmap(
		0:M_dim, 0:M_dim, P_mat,
		yflip = true,
		xticks = (0:M_dim, etiquetas_estados), yticks = (0:M_dim, etiquetas_estados),
		title = "Matriz de transición P — $(M_dim + 1)×$(M_dim + 1) (λ = $(@sprintf("%.2f", lambda_p)))",
		xlabel = "Estado destino en la siguiente generación (j)",
		ylabel = "Estado origen actual (i)",
		colorbar_title = "Probabilidad",
		color = cmap_mat_markov, clims = (0.0, 1.0),
		size = (680, 560), aspect_ratio = :equal
	)

	# Valor de cada entrada P(i,j), centrado en su celda
	for i in 0:M_dim
		for j in 0:M_dim
			txt = @sprintf("%.2f", P_mat[i + 1, j + 1])
			annotate!(p_heatmap_P, j, i, text(txt, font(fsize_p, :center, :black, fuente_negrita)))
		end
	end

	p_heatmap_P
end

# ╔═╡ 00000002-0002-0002-0002-000000000065
begin
	# Verificación de que P es una matriz estocástica: cada fila suma 1
	error_filas = maximum(abs.(vec(sum(P_mat, dims = 2)) .- 1.0))
	error_txt = @sprintf("%.1e", error_filas)

	md"""
	##### Interpretación directa de la matriz

	- **Fila 0.** Es $[1.00, 0.00, 0.00, \ldots]$ y representa el **pozo absorbente**: una vez que el apellido llega a 0 portadores, nunca sale.
	- **Filas $i > 0$.** Entrada $\left(i,j\right)$ = probabilidad de pasar de $i$ a $j$ personas en la siguiente generación. Cada fila suma 1 (error numérico máximo: $(error_txt)).
	- **Última fila y última columna ("$M+$").** La última columna acumula la probabilidad de alcanzar $M$ o más portadores. Los valores menores que 0.005 se muestran como 0.00.
	- **Sesgo hacia la izquierda.** Cuando $\lambda < 1$, la masa de cada fila se concentra a la izquierda de la diagonal (hacia números menores de portadores), lo que empuja sistemáticamente la distribución hacia el estado 0. Si $\lambda > 1$, la masa se desplaza hacia la derecha.
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000067
begin
	M_curr = M_slider
	z0_c = min(z0_markov_slider, M_curr)     # el estado inicial no puede superar el tope M

	# Distribución inicial μ₀ = δ_{z0} (sección 2.2): vector fila con un 1 en la
	# posición z0 (índice z0+1) y ceros en las demás
	pi_0 = zeros(Float64, M_curr + 1)
	pi_0[z0_c + 1] = 1.0

	# Distribución en la generación n: πₙ = π₀ Pⁿ (sección 2.5). Se aplica P una vez
	# por generación: πₖ₊₁ = πₖ P (vector fila por matriz).
	pi_n = copy(pi_0)
	for _ in 1:gen_pi_slider
		pi_n = vec(pi_n' * P_mat)
	end

	colores_pi = [k == 0 ? "#e53935" : "#1e88e5" for k in 0:M_curr]

	# Ejes numéricos con marcas de texto: las anotaciones quedan centradas sobre las barras
	p_bar_pi = bar(
		0:M_curr, pi_n,
		xticks = (0:M_curr, etiquetas_estados),
		xlabel = "Número de portadores (estado)",
		ylabel = "Probabilidad πₙ(k)",
		title = "Distribución πₙ en la generación n = $(gen_pi_slider) (Z₀ = $(z0_c))",
		color = colores_pi, alpha = 0.85, legend = false,
		ylims = (0.0, 1.15), size = (800, 400)
	)

	# Porcentaje sobre las barras con masa apreciable
	for (k, p_val) in zip(0:M_curr, pi_n)
		if p_val > 0.015
			annotate!(p_bar_pi, k, p_val + 0.04, text(@sprintf("%.1f%%", p_val * 100), 8, :center, :black, fuente_negrita))
		end
	end

	p_bar_pi
end

# ╔═╡ 00000002-0002-0002-0002-000000000068
begin
	max_g_flujo = 20
	M_flujo = M_slider
	z0_f = min(z0_markov_slider, M_flujo)

	# Cada fila de mat_flujo es πₙ para n = 0, 1, …, max_g_flujo (πₙ₊₁ = πₙ P)
	mat_flujo = zeros(Float64, max_g_flujo + 1, M_flujo + 1)
	vec_temp = zeros(Float64, M_flujo + 1)
	vec_temp[z0_f + 1] = 1.0
	mat_flujo[1, :] = vec_temp

	for g in 1:max_g_flujo
		vec_temp = vec(vec_temp' * P_mat)
		mat_flujo[g + 1, :] = vec_temp
	end

	p_flujo = heatmap(
		0:M_flujo, 0:max_g_flujo, mat_flujo,
		yflip = true,
		xticks = (0:M_flujo, etiquetas_estados), yticks = 0:2:max_g_flujo,
		title = "Evolución del flujo de probabilidad desde Z₀ = $(z0_f)",
		xlabel = "Estado (portadores)", ylabel = "Generación (n)",
		colorbar_title = "Probabilidad",
		color = cmap_mat_markov, clims = (0.0, 1.0),
		size = (760, 540)
	)
	p_flujo
end

# ╔═╡ 00000002-0002-0002-0002-000000000073
begin
	# Modelo de juguete para el grafo: estados 0, 1, 2, 3+ (M = 3), con la misma
	# función que construye la matriz de la sección anterior
	M_toy = 3
	P_toy = matriz_transicion_truncada(lambda_col_slider, M_toy)

	# ── Geometría, en coordenadas de datos (aspect_ratio = :equal da círculos reales) ──
	xs_toy = [0.0, 3.6, 7.2, 10.8]                  # centros de los nodos (todos en y = 0)
	r_nodo = 0.62                                   # radio de los nodos
	nombres_toy = ["0", "1", "2", "3+"]
	colores_nodos = ["#e53935", "#1e88e5", "#1e88e5", "#8e24aa"]
	umbral_prob = 0.025                             # no se dibujan transiciones menos probables

	# Alturas de los arcos según la distancia d = |j - i| entre nodos
	altura_arco_sup(d) = d == 1 ? 0.95 : 2.45 + 0.6 * (d - 2)   # los de d ≥ 2 pasan por encima de los autolazos
	altura_arco_inf(d) = 0.95 + 0.6 * (d - 1)
	# Ángulo (en grados) con que el arco sale/llega al borde del nodo; varía con d para
	# que las puntas de flecha que llegan al mismo nodo no se superpongan. Los arcos
	# inferiores usan un paso mayor porque en esa mitad no hay autolazos.
	angulo_borde_sup(d) = 22.0 + 22.0 * (d - 1)
	angulo_borde_inf(d) = 22.0 + 22.0 * (d - 1)

	# Cabeza de flecha triangular con la punta en (tx, ty) y dirección unitaria (ux, uy)
	function cabeza_flecha(tx, ty, ux, uy; largo = 0.26, semiancho = 0.12)
		bx, by = tx - ux * largo, ty - uy * largo
		return Shape([tx, bx - uy * semiancho, bx + uy * semiancho],
		             [ty, by + ux * semiancho, by - ux * semiancho])
	end

	# Arco de Bézier cuadrático entre dos puntos del borde de los nodos, con flecha
	# en el destino. `h` es la altura (con signo) del vértice del arco.
	function dibujar_arco!(plt, xs, ys, xe, ye, h, col, lw, prob, etiqueta_dentro)
		cx = (xs + xe) / 2
		cy = 2 * h - ys                              # punto de control: el vértice del arco queda en y = h
		bez(t) = ((1 - t)^2 * xs + 2 * (1 - t) * t * cx + t^2 * xe,
		          (1 - t)^2 * ys + 2 * (1 - t) * t * cy + t^2 * ye)
		largo_flecha = 0.26
		# Se recorta la curva donde empieza la cabeza de la flecha
		t_corte = 1.0
		for t in range(1.0, 0.5, length = 200)
			px, py = bez(t)
			if hypot(px - xe, py - ye) >= 0.8 * largo_flecha
				t_corte = t
				break
			end
		end
		ts_c = range(0.0, t_corte, length = 60)
		pts = [bez(t) for t in ts_c]
		plot!(plt, first.(pts), last.(pts), color = col, linewidth = lw, alpha = 0.85, label = false)
		# La tangente en t = 1 es (E − control)
		tx, ty = xe - cx, ye - cy
		L = hypot(tx, ty)
		plot!(plt, cabeza_flecha(xe, ye, tx / L, ty / L; largo = largo_flecha),
		      color = col, linecolor = col, linewidth = 1, alpha = 0.95, label = false)
		# Etiqueta en el vértice del arco, del lado de dentro (hacia el eje) o de fuera
		signo = h > 0 ? 1.0 : -1.0
		desplaza = etiqueta_dentro ? -0.27 * signo : 0.30 * signo
		annotate!(plt, cx, h + desplaza,
		          text(@sprintf("%.2f", prob), 8, :center, col, fuente_negrita))
	end

	# Autolazo sobre el nodo (cx, cy): circunferencia de radio r_l que corta el borde
	# del nodo, recorrida en sentido antihorario, con la flecha entrando al nodo
	function dibujar_autolazo!(plt, cx, cy, col, lw, prob)
		d_l, r_l = 0.80, 0.42
		# Ángulo θ0 donde la circunferencia del lazo corta el borde del nodo:
		# |(r_l cos θ, d_l + r_l sin θ)| = r_nodo  ⇒  sin θ0 = (r_nodo² − r_l² − d_l²) / (2 d_l r_l)
		θ0 = asin((r_nodo^2 - r_l^2 - d_l^2) / (2 * d_l * r_l))
		θ1 = π - θ0                                   # extremo final (lado izquierdo del nodo)
		largo_flecha = 0.22
		θ_corte = θ1 - 0.8 * largo_flecha / r_l       # recorte de la curva donde empieza la flecha
		thetas = range(θ0, θ_corte, length = 60)
		plot!(plt, cx .+ r_l .* cos.(thetas), cy .+ d_l .+ r_l .* sin.(thetas),
		      color = col, linewidth = lw, alpha = 0.85, label = false)
		# Punta de flecha en el extremo θ1, tangente al sentido antihorario (−sin θ, cos θ)
		tipx, tipy = cx + r_l * cos(θ1), cy + d_l + r_l * sin(θ1)
		plot!(plt, cabeza_flecha(tipx, tipy, -sin(θ1), cos(θ1); largo = largo_flecha, semiancho = 0.11),
		      color = col, linecolor = col, linewidth = 1, alpha = 0.95, label = false)
		annotate!(plt, cx, cy + d_l + r_l + 0.30,
		          text(@sprintf("%.2f", prob), 8, :center, col, fuente_negrita))
	end

	# ── Selección de las transiciones a dibujar y límites verticales de la figura ──
	hay_arco_sup_largo = any(P_toy[i + 1, j + 1] > umbral_prob for i in 0:M_toy for j in 0:M_toy if j - i >= 2)
	d_inf_max = maximum([i - j for i in 0:M_toy for j in 0:M_toy if i > j && P_toy[i + 1, j + 1] > umbral_prob]; init = 1)
	y_top = hay_arco_sup_largo ? 3.2 : 2.3
	y_bot = -(altura_arco_inf(d_inf_max) + 0.55)
	x_lo, x_hi = -1.5, 12.3

	p_diag_markov = plot(
		xlims = (x_lo, x_hi), ylims = (y_bot, y_top),
		aspect_ratio = :equal,
		title = "Grafo dirigido de transiciones (M = 3, λ = $(@sprintf("%.2f", lambda_col_slider)))",
		legend = false, grid = false, framestyle = :none,
		size = (940, round(Int, 940 * (y_top - y_bot) / (x_hi - x_lo)) + 70)
	)

	# Se dibujan primero las transiciones y luego los nodos (que tapan los extremos)
	for i in 0:M_toy
		for j in 0:M_toy
			prob = P_toy[i + 1, j + 1]
			prob > umbral_prob || continue
			if i == j
				# Autolazo: rojo grueso en el estado absorbente, gris en los demás
				dibujar_autolazo!(p_diag_markov, xs_toy[i + 1], 0.0,
				                  i == 0 ? "#e53935" : "#757575", i == 0 ? 2.5 : 1.5, prob)
			elseif j > i
				# Crecimiento: arco superior azul
				d = j - i
				a = deg2rad(angulo_borde_sup(d))
				xs, ys = xs_toy[i + 1] + r_nodo * cos(a), r_nodo * sin(a)
				xe, ye = xs_toy[j + 1] + r_nodo * cos(π - a), r_nodo * sin(π - a)
				dibujar_arco!(p_diag_markov, xs, ys, xe, ye, altura_arco_sup(d),
				              "#1e88e5", 1.2 + 3.0 * prob, prob, d == 1)
			else
				# Decaimiento: arco inferior; rojo si llega a la extinción, naranja en otro caso
				d = i - j
				a = deg2rad(angulo_borde_inf(d))
				xs, ys = xs_toy[i + 1] + r_nodo * cos(π + a), -r_nodo * sin(a)
				xe, ye = xs_toy[j + 1] + r_nodo * cos(-a), -r_nodo * sin(a)
				dibujar_arco!(p_diag_markov, xs, ys, xe, ye, -altura_arco_inf(d),
				              j == 0 ? "#e53935" : "#e65100", 1.2 + 3.5 * prob, prob, true)
			end
		end
	end

	# Nodos: círculos blancos con borde de color y etiqueta centrada
	θc = range(0, 2π, length = 80)
	for idx in 1:(M_toy + 1)
		x0 = xs_toy[idx]
		plot!(p_diag_markov, Shape(x0 .+ r_nodo .* cos.(θc), r_nodo .* sin.(θc)),
		      fillcolor = :white, linecolor = colores_nodos[idx], linewidth = 3.5, label = false)
		if idx == 1
			annotate!(p_diag_markov, x0, 0.16, text("0", 13, :center, :black, fuente_negrita))
			annotate!(p_diag_markov, x0, -0.24, text("extinto", 8, :center, :black, fuente_negrita))
		else
			annotate!(p_diag_markov, x0, 0.0, text(nombres_toy[idx], 13, :center, :black, fuente_negrita))
		end
	end

	p_diag_markov
end

# ╔═╡ 00000002-0002-0002-0002-000000000043
md"**Evaluar la probabilidad de extinción acumulada hasta la generación ($X$):** $(@bind gen_eval_slider Slider(5:1:40, default=20, show_value=true))"

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

	# Población inicial Z₀ = número de portadores del apellido seleccionado
	z0_seleccionado = dane_dict[ap_col_select]

	# Modelo de descendencia: Poisson(λ) para los hijos varones de cada portador,
	# instanciado con ModeloPBGW de src/markovModels.jl
	modelo_colombia = ModeloPBGW(lambda_col_slider; max_hijos = 15)

	# Tabla comparativa que resalta el apellido activo
	filas_tabla = [
		"| Puesto | Apellido | Portadores DANE 2018 (\$Z_0\$) | Estado en la simulación |",
		"|:---|:---|:---|:---|"
	]
	for (idx, (ap, cnt)) in enumerate(zip(apellidos_colombia, conteos_dane))
		marca = ap == ap_col_select ? "**SELECCIONADO**" : "En espera"
		ap_fmt = ap == ap_col_select ? "**$ap**" : ap
		push!(filas_tabla, "| $idx | $ap_fmt | $(fmt_miles(cnt)) | $marca |")
	end
	Markdown.parse(join(filas_tabla, "\n"))
end

# ╔═╡ 00000002-0002-0002-0002-000000000069
begin
	gens_curva_p = 0:25
	z0_curva = min(z0_markov_slider, M_dim)

	# Curva de la cadena truncada: P(Zₙ = 0) es la componente 0 de πₙ = π₀ Pⁿ
	prob_ext_cadena = Float64[]
	vec_t = zeros(Float64, M_dim + 1)
	vec_t[z0_curva + 1] = 1.0
	for g in gens_curva_p
		push!(prob_ext_cadena, vec_t[1])          # P(Z_g = 0)
		vec_t = vec(vec_t' * P_mat)
	end

	# Curva exacta del proceso sin truncar: P(Zₙ = 0) = [gₙ(0)]^Z₀ (sección 4.3)
	prob_ext_exacta = prob_extincion_curva(modelo_colombia, last(gens_curva_p), z0_curva)

	# La leyenda va debajo del área de datos para no tapar ninguna de las curvas
	p_ext_acum = plot(
		gens_curva_p, prob_ext_cadena,
		xlabel = "Generación (n)", ylabel = "Probabilidad acumulada P(Zₙ = 0)",
		title = "Curva de absorción: ¿con qué rapidez cae la cadena al estado 0?",
		label = "Cadena truncada (matriz P)", color = :red, linewidth = 3,
		size = (800, 460), ylims = (0.0, 1.05), grid = true,
		legend = :outerbottom, legend_columns = 2
	)
	plot!(p_ext_acum, gens_curva_p, prob_ext_exacta,
		label = "Proceso exacto: [gₙ(0)]^Z₀", color = :black, linewidth = 2, linestyle = :dash)
	hline!(p_ext_acum, [0.5], label = "Mediana (50 %)", color = :gray, linestyle = :dash)
	hline!(p_ext_acum, [1.0], label = "Absorción total (100 %)", color = :black, linestyle = :dot)
	p_ext_acum
end

# ╔═╡ 00000002-0002-0002-0002-000000000047
md"""
### Aspectos numéricos y tipo de Monte Carlo utilizado

#### 1. Cómo se simula y qué aporta la propiedad reproductiva de Poisson

En `src/markovModels.jl`, `simular_generaciones` recorre los individuos uno por uno y sortea el número de hijos de cada uno con la distribución categórica `Categorical(modelo.pk)`:

```julia
for _ in 1:poblacion
    hijos_totales += rand(rng, dist_hijos) - 1
end
```

Si la descendencia individual es Poisson($\lambda$), la **propiedad reproductiva** (la suma de Poisson independientes es Poisson con la suma de las medias) implica

$$\sum_{i=1}^{Z_n} X_i \sim \mathrm{Poisson}(Z_n \lambda).$$

Esta identidad tiene dos usos:

* **Construcción exacta de $P$.** La fila $i$ de la matriz de transición es una Poisson$\left(i\lambda\right)$ (sección 7.2), por lo que la matriz de la Parte IV se obtiene sin convoluciones.
* **Simulación más rápida (no implementada).** Permitiría generar cada generación con un único sorteo Poisson$\left(Z_n\lambda\right)$, de costo constante. La implementación actual suma individuo por individuo: es exacta y equivalente en distribución, pero su costo crece con el tamaño de la población, del orden de $10^5$ a $10^6$ portadores en esta sección. Por eso los cálculos de esta parte son más lentos que los de la Parte II.
"""

# ╔═╡ 00000002-0002-0002-0002-000000000048
begin
	# Simulación Monte Carlo con las funciones de src/markovModels.jl.
	# Devuelve las estadísticas del tiempo de extinción y una muestra de trayectorias.
	function correr_montecarlo_colombia(modelo::ModeloPBGW, z0::Int, n_sims::Int, gen_eval::Int, semilla::Int)
		rng = MersenneTwister(semilla)
		# Horizonte de simulación: al menos 35 generaciones y 15 más allá de la evaluada
		max_gen_sim = max(gen_eval + 15, 35)
		tiempos_ext = Int[]            # generación de extinción de cada simulación que se extinguió
		trayectorias_muestra = Vector{Vector{Int}}()
		extintos_en_X = 0              # simulaciones extintas en ≤ gen_eval generaciones

		for s in 1:n_sims
			traj = simular_generaciones(modelo, max_gen_sim; poblacion_inicial = z0, rng = rng)
			if s <= 6
				push!(trayectorias_muestra, copy(traj))
			end
			# traj[k] es la población en la generación k-1; el primer 0 marca la extinción
			idx_cero = findfirst(==(0), traj)
			if idx_cero !== nothing
				gen_ext = idx_cero - 1
				push!(tiempos_ext, gen_ext)
				if gen_ext <= gen_eval
					extintos_en_X += 1
				end
			end
		end

		# Estimador de P(T ≤ gen_eval): fracción de simulaciones extintas hasta gen_eval
		prob_ext_X = extintos_en_X / n_sims

		# Simulaciones censuradas: no se extinguieron dentro del horizonte simulado.
		# Los estadísticos del tiempo de extinción solo son válidos si se tiene en cuenta
		# esta censura, porque promediar únicamente las simulaciones extintas sesga a la baja.
		n_censuradas = n_sims - length(tiempos_ext)

		# Mediana de T entre las n_sims simulaciones: el elemento k_med-ésimo en orden
		# creciente (las censuradas ocupan las últimas posiciones). Solo está determinada
		# si al menos la mitad de las simulaciones se extinguió.
		ordenados = sort(tiempos_ext)
		k_med = cld(n_sims, 2)
		mediana_ext = length(ordenados) >= k_med ? ordenados[k_med] : nothing

		# La media solo está determinada si ninguna simulación quedó censurada
		media_ext = n_censuradas == 0 ? sum(tiempos_ext) / n_sims : nothing

		return (
			tiempos_ext = tiempos_ext,
			trayectorias_muestra = trayectorias_muestra,
			extintos_en_X = extintos_en_X,
			prob_ext_X = prob_ext_X,
			n_censuradas = n_censuradas,
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
	n_censuradas = res_sim_col.n_censuradas
	mediana_ext = res_sim_col.mediana_ext
	media_ext = res_sim_col.media_ext
	max_gen_sim = res_sim_col.max_gen_sim

	# ── Validación teórica con las funciones del proyecto ──
	m_esp = esperanza(modelo_colombia)                    # m = λ (salvo truncamiento)
	q_teorico = probabilidad_extincion(modelo_colombia)   # q = 1 si m ≤ 1

	# Aproximación determinista (término dominante de la sección 4.4):
	# E[Z_n] = Z₀ mⁿ = 1  ⇒  n ≈ ln(Z₀) / ln(1/m)
	t_teorico = m_esp < 1.0 ? log(z0_seleccionado) / log(1.0 / m_esp) : Inf

	# Teoría exacta (sección 4.3): P(T ≤ n) = [g_n(0)]^Z₀, evaluada por iteración de fgp
	prob_ext_X_teo = prob_extincion_hasta(modelo_colombia, gen_eval_slider, z0_seleccionado)
	mediana_teo = mediana_extincion_teorica(modelo_colombia, z0_seleccionado)

	# Textos de las filas de la tabla (se construyen aparte para tratar los casos sin valor)
	txt_mediana_mc = mediana_ext === nothing ?
		"No observada: más de la mitad de las simulaciones sigue viva tras $(max_gen_sim) generaciones" :
		"Gen $(mediana_ext) (~$(mediana_ext * anios_por_generacion) años)"
	txt_mediana_teo = mediana_teo === nothing ?
		"No alcanzada en el rango calculado" :
		"Gen $(mediana_teo) (~$(mediana_teo * anios_por_generacion) años)"
	txt_media_mc = media_ext === nothing ?
		"No determinada: $(n_censuradas) simulaciones censuradas" :
		"Gen $(round(media_ext, digits=1))"
	txt_regimen = clasificar_regimen(m_esp)

	md"""
	### Resultados del Monte Carlo para **$(ap_col_select)**

	| Métrica | Valor | Significado y validación teórica |
	|:---|:---|:---|
	| **Población inicial (Z₀)** | **$(fmt_miles(z0_seleccionado))** | Portadores del apellido según la tabla de datos |
	| **Media de descendencia (m = E[X])** | **$(round(m_esp, digits=3))** | Evaluada con `esperanza(modelo)` |
	| **Régimen del proceso** | **$(txt_regimen)** | Evaluado con `es_extincion_segura(modelo)` |
	| **Probabilidad asintótica de extinción (q)** | **$(round(q_teorico * 100, digits=2))%** | Calculada con `probabilidad_extincion(modelo)` |
	| **Extinción hasta la generación $(gen_eval_slider), Monte Carlo** | **$(round(prob_ext_X * 100, digits=2))%** | Fracción observada con `simular_generaciones` |
	| **Extinción hasta la generación $(gen_eval_slider), teoría exacta** | **$(round(prob_ext_X_teo * 100, digits=2))%** | [gₙ(0)]^Z₀ por iteración de `fgp` (sección 4.3) |
	| **Generación mediana de extinción, Monte Carlo** | **$(txt_mediana_mc)** | El 50 % de las simulaciones se extinguió en esa generación o antes |
	| **Generación mediana de extinción, teoría exacta** | **$(txt_mediana_teo)** | Menor n con [gₙ(0)]^Z₀ ≥ 1/2 |
	| **Generación media de extinción, Monte Carlo** | **$(txt_media_mc)** | Promedio aritmético de los tiempos de extinción |
	| **Aproximación determinista** | **Gen $(round(t_teorico, digits=1))** | Término dominante T ≈ ln(Z₀) / ln(1/m) de la sección 4.4 |
	"""
end

# ╔═╡ 00000002-0002-0002-0002-000000000049
begin
	# Distribución del tiempo de extinción: frecuencias observadas frente a la
	# distribución teórica exacta
	if !isempty(tiempos_ext)
		g_min, g_max = minimum(tiempos_ext), maximum(tiempos_ext)
		gens_x = g_min:g_max
		# Número de simulaciones que se extinguieron exactamente en la generación g
		conteos = [count(==(g), tiempos_ext) for g in gens_x]

		# Teoría: n_sims · P(T = g) = n_sims · [P(T ≤ g) − P(T ≤ g−1)], con
		# P(T ≤ g) = [g_g(0)]^Z₀  (la posición g+1 de la curva corresponde a la generación g)
		curva_acum = prob_extincion_curva(modelo_colombia, g_max, z0_seleccionado)
		teorico = [n_sims_col_slider * (curva_acum[g + 1] - (g >= 1 ? curva_acum[g] : 0.0)) for g in gens_x]

		# Las barras hasta la generación X se resaltan en rojo; el resto en gris
		colores_barras_hist = [g <= gen_eval_slider ? "#e53935" : "#b0bec5" for g in gens_x]
		paso_x = max(1, cld(length(gens_x), 12))           # marcas enteras y legibles
		y_tope = 1.3 * max(maximum(conteos), maximum(teorico))

		p_hist_ext = bar(
			collect(gens_x), conteos,
			bar_width = 1.0, color = colores_barras_hist, alpha = 0.85, label = false,
			xlabel = "Generación de extinción (1 generación ≈ $(anios_por_generacion) años)",
			ylabel = "Número de simulaciones",
			title = "Tiempo de extinción de $(ap_col_select) (N = $(n_sims_col_slider))",
			legend = :topright, size = (900, 520),
			xticks = collect(g_min:paso_x:g_max), ylims = (0, y_tope),
			xlims = (g_min - 1.5, g_max + 1.5)
		)
		plot!(p_hist_ext, collect(gens_x), teorico,
			label = "Teoría exacta: N · P(T = k)",
			color = :black, marker = :circle, markersize = 4, linewidth = 1.5)
		# Línea de la mediana observada (si está determinada)
		if mediana_ext !== nothing
			vline!(p_hist_ext, [mediana_ext], label = "Mediana observada = Gen $(mediana_ext)",
				color = :black, linewidth = 2, linestyle = :dash)
		end
		# El umbral "≤ X" corresponde al borde derecho de la barra X, es decir, X + 0.5
		vline!(p_hist_ext, [gen_eval_slider + 0.5],
			label = "Umbral X = Gen $(gen_eval_slider): $(round(prob_ext_X * 100, digits=1)) % extinto",
			color = :blue, linewidth = 2)
		p_hist_ext
	else
		# Ninguna simulación se extinguió dentro del horizonte (censura total)
		txt_teo = mediana_teo === nothing ? "no alcanzada en el rango calculado" : "Gen $(mediana_teo)"
		plot(title = "Ninguna simulación se extinguió en $(max_gen_sim) generaciones (mediana teórica exacta: $(txt_teo))",
			titlefontsize = 11, framestyle = :none, size = (900, 200))
	end
end

# ╔═╡ 00000002-0002-0002-0002-000000000050
begin
	# Trayectorias en escala logarítmica. En esta escala Zₙ = 0 no es representable;
	# las extinciones se dibujan en piso_log = 0.5.
	piso_log = 0.5
	gens_x_log = 0:max_gen_sim

	# Media esperada E[Zₙ] = Z₀ λⁿ (sección 4.2), sin recortes: en escala log es una recta
	curva_esperada = [z0_seleccionado * lambda_col_slider^g for g in gens_x_log]
	y_max_log = 3 * max(z0_seleccionado, maximum(maximum.(trayectorias_muestra)))
	# La leyenda va en la esquina libre: arriba a la derecha si el linaje decae
	pos_leyenda_log = lambda_col_slider < 1.0 ? :topright : :bottomright

	# tope_poblacion por defecto de simular_generaciones (src/markovModels.jl): un linaje
	# que lo alcanza se declara establecido y su población se congela en ese valor
	tope_pob_col = 1_000_000
	tit_log = lambda_col_slider < 1.0 ? "Decaimiento estocástico de linajes" : "Trayectorias estocásticas de linajes"

	p_log = plot(
		xlabel = "Generación (n)", ylabel = "Portadores (escala log₁₀)",
		title = "$(tit_log) — $(ap_col_select) ($(@sprintf("%.2f", lambda_col_slider)) varones por varón)",
		yaxis = :log10, size = (850, 440), legend = pos_leyenda_log,
		ylims = (0.3, y_max_log),
		grid = true, gridstyle = :dash, gridalpha = 0.3
	)

	plot!(p_log, gens_x_log, curva_esperada, label = "Media esperada E[Zₙ] = Z₀·λⁿ",
		color = :black, linewidth = 2.5, linestyle = :dash)

	for (idx, tr) in enumerate(trayectorias_muestra)
		tr_plot = [max(piso_log, Float64(x)) for x in tr]     # los ceros se dibujan en piso_log
		plot!(p_log, 0:(length(tr_plot) - 1), tr_plot,
			label = idx == 1 ? "Simulaciones individuales (Zₙ = 0 se dibuja en 0.5)" : false,
			alpha = 0.7, linewidth = 1.5)
	end
	hline!(p_log, [1.0], label = "1 portador", color = :red, linestyle = :dot, linewidth = 1.5)
	if any(maximum(tr) >= tope_pob_col for tr in trayectorias_muestra)
		hline!(p_log, [tope_pob_col], label = "Tope de simulación (10⁶): linaje establecido",
			color = :gray, linestyle = :dot, linewidth = 1.5)
	end
	p_log
end

# ╔═╡ 00000002-0002-0002-0002-000000000051
begin
	# Comparación de los 10 apellidos: mediana del tiempo de extinción obtenida con
	# simular_generaciones (n_rep simulaciones por apellido, horizonte gen_max).
	function calcular_medianas_10(modelo::ModeloPBGW; gen_max::Int = 35, n_rep::Int = 30)
		medianas = Float64[]
		rng_comp = MersenneTwister(123)
		for cnt in conteos_dane
			tiempos_i = Int[]
			for _ in 1:n_rep
				traj_i = simular_generaciones(modelo, gen_max; poblacion_inicial = cnt, rng = rng_comp)
				idx_0 = findfirst(==(0), traj_i)
				# Una simulación que no se extingue dentro del horizonte se registra como
				# gen_max + 1 (censura): la mediana resultante es solo una cota inferior.
				push!(tiempos_i, idx_0 !== nothing ? idx_0 - 1 : gen_max + 1)
			end
			push!(medianas, median(tiempos_i))
		end
		return medianas
	end

	gen_max_comp = 35
	tiempos_medios_10 = calcular_medianas_10(modelo_colombia; gen_max = gen_max_comp)
	censurado_10 = tiempos_medios_10 .> gen_max_comp          # mediana fuera del horizonte simulado

	# Referencia teórica exacta para cada apellido (iteración de la FGP, sección 4.3)
	medianas_teo_10 = [mediana_extincion_teorica(modelo_colombia, cnt) for cnt in conteos_dane]

	# Colores: apellido seleccionado en rojo, el resto en azul, censurados en gris
	colores_barras = [censurado_10[i] ? "#9e9e9e" : (apellidos_colombia[i] == ap_col_select ? "#e53935" : "#1e88e5") for i in 1:10]

	# Orden inverso para que el puesto 1 (Rodríguez) aparezca arriba
	posiciones = 1:10
	aps_rev = reverse(apellidos_colombia)
	tiempos_rev = reverse(tiempos_medios_10)
	censurado_rev = reverse(censurado_10)
	colores_rev = reverse(colores_barras)
	teo_rev = reverse(medianas_teo_10)

	# Marcadores de la mediana teórica exacta (donde está definida)
	pos_teo = Int[posiciones[i] for i in 1:10 if teo_rev[i] !== nothing]
	val_teo = Float64[teo_rev[i] for i in 1:10 if teo_rev[i] !== nothing]

	# El eje X debe abarcar tanto las barras como las medianas teóricas
	x_tope = maximum([tiempos_medios_10; val_teo]) * 1.45

	p_comp = bar(
		posiciones, tiempos_rev,
		orientation = :horizontal,
		yticks = (posiciones, aps_rev),
		xlabel = "Generaciones hasta la extinción (mediana de 30 simulaciones)",
		ylabel = "Apellido",
		title = "¿Cuánto tarda en extinguirse cada apellido? (λ = $(@sprintf("%.2f", lambda_col_slider)))",
		color = colores_rev, alpha = 0.85, label = false,
		legend = isempty(val_teo) ? false : :outerbottom, size = (850, 500),
		ylims = (0.4, 10.6), xlims = (0, x_tope)
	)
	scatter!(p_comp, val_teo, pos_teo, marker = :diamond, markersize = 6, color = :black,
		label = "Mediana teórica exacta (iteración de la FGP)")

	# Etiqueta de cada barra, a la derecha del extremo más lejano (barra o marcador)
	for (idx, t_med) in enumerate(tiempos_rev)
		x_fin = teo_rev[idx] === nothing ? t_med : max(t_med, teo_rev[idx])
		etiqueta = if censurado_rev[idx]
			"> $(gen_max_comp) gen. (censurado)"
		else
			"Gen $(round(Int, t_med)) (~$(round(Int, t_med * anios_por_generacion)) años)"
		end
		annotate!(p_comp, x_fin + 0.025 * x_tope, idx, text(etiqueta, 8, :left, :black, fuente_negrita))
	end
	p_comp
end

# ╔═╡ 00000002-0002-0002-0002-000000000052
begin
	# Comparación entre el apellido más numeroso y el menos numeroso de la lista
	i_max = argmax(conteos_dane)
	i_min = argmin(conteos_dane)
	z_mayor, z_menor = conteos_dane[i_max], conteos_dane[i_min]
	ap_mayor, ap_menor = apellidos_colombia[i_max], apellidos_colombia[i_min]
	lam_txt = @sprintf("%.2f", lambda_col_slider)

	if lambda_col_slider < 1.0
		# Mediana exacta de extinción de cada uno (iteración de la FGP)
		med_mayor = mediana_extincion_teorica(modelo_colombia, z_mayor)
		med_menor = mediana_extincion_teorica(modelo_colombia, z_menor)
		# Generaciones adicionales que aporta multiplicar la población: ln(k) / ln(1/λ)
		extra_gen = log(z_mayor / z_menor) / log(1 / lambda_col_slider)
		extra_txt = @sprintf("%.2f", extra_gen)
		# Textos precalculados: dentro de md solo se interpolan variables simples
		z_mayor_txt = fmt_miles(z_mayor)
		z_menor_txt = fmt_miles(z_menor)
		dif_txt = fmt_miles(z_mayor - z_menor)
		razon_txt = string(round(z_mayor / z_menor, digits = 2))
		dif_gen = (med_mayor === nothing || med_menor === nothing) ? "?" : string(med_mayor - med_menor)  # diferencia entre medianas

		# Reglas de interpolación en md: sin paréntesis literal justo antes de $( y sin
		# negritas alrededor de una interpolación que siga a otra en el mismo párrafo
		md"""
		!!! tip "La paradoja de los apellidos masivos"
		    Con λ = $(lam_txt), la mediana exacta de extinción es de $(med_mayor) generaciones para $(ap_mayor), con $(z_mayor_txt) personas, y de $(med_menor) generaciones para $(ap_menor), con $(z_menor_txt) personas. Hay $(dif_txt) portadores de diferencia entre ambos y, sin embargo, la diferencia en la mediana del tiempo de extinción es de solo $(dif_gen) generación(es).

		    ¿Por qué? Porque en un proceso subcrítico la población decae **geométricamente**, E[Zₙ] = Z₀·λⁿ. Al despejar el tiempo en que la población esperada llega a un portador:

		    $$T_{\text{ext}} \approx \frac{\ln Z_0}{\ln(1/\lambda)}$$

		    El tiempo crece con el **logaritmo** de la población inicial: multiplicar los portadores por un factor k añade solo ln(k) / ln(1/λ) generaciones. Aquí, el factor $(razon_txt) entre ambos apellidos equivale a $(extra_txt) generaciones adicionales.
		"""
	else
		md"""
		!!! note "Dependencia de la población inicial"
		    Con λ = $(lam_txt) ≥ 1 el proceso no es subcrítico y la fórmula logarítmica del tiempo de extinción no aplica. En el caso crítico (λ = 1) la extinción es segura pero muy lenta; en el supercrítico (λ > 1) la extinción tiene probabilidad menor que 1 y, para poblaciones de cientos de miles de portadores, es prácticamente imposible (sección 5).
		"""
	end
end

# ╔═╡ 4c85c2e9-5396-43d7-bea1-aed1414de316
html"""
<style>
    .markdown p, .markdown li {
        text-align: justify;
    }
</style>
"""

# ╔═╡ Cell order:
# ╟─00000002-0002-0002-0002-000000000001
# ╟─00000002-0002-0002-0002-000000000085
# ╠═00000002-0002-0002-0002-000000000002
# ╟─00000002-0002-0002-0002-000000000074
# ╟─00000002-0002-0002-0002-000000000075
# ╟─00000002-0002-0002-0002-000000000084
# ╟─00000002-0002-0002-0002-000000000062
# ╟─00000002-0002-0002-0002-000000000063
# ╟─00000002-0002-0002-0002-000000000083
# ╠═00000002-0002-0002-0002-000000000064
# ╟─00000002-0002-0002-0002-000000000065
# ╟─00000002-0002-0002-0002-000000000072
# ╟─00000002-0002-0002-0002-000000000073
# ╟─00000002-0002-0002-0002-000000000066
# ╟─00000002-0002-0002-0002-000000000067
# ╟─00000002-0002-0002-0002-000000000068
# ╟─00000002-0002-0002-0002-000000000069
# ╠═00000002-0002-0002-0002-000000000076
# ╟─00000002-0002-0002-0002-000000000077
# ╟─00000002-0002-0002-0002-000000000003
# ╟─00000002-0002-0002-0002-000000000004
# ╟─00000002-0002-0002-0002-000000000078
# ╟─00000002-0002-0002-0002-000000000079
# ╟─00000002-0002-0002-0002-000000000080
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
# ╠═00000002-0002-0002-0002-000000000022
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
# ╠═00000002-0002-0002-0002-000000000042
# ╟─00000002-0002-0002-0002-000000000043
# ╟─00000002-0002-0002-0002-000000000044
# ╟─00000002-0002-0002-0002-000000000045
# ╟─00000002-0002-0002-0002-000000000046
# ╠═00000002-0002-0002-0002-000000000047
# ╟─00000002-0002-0002-0002-000000000048
# ╟─00000002-0002-0002-0002-000000000049
# ╟─00000002-0002-0002-0002-000000000050
# ╟─00000002-0002-0002-0002-000000000051
# ╟─00000002-0002-0002-0002-000000000052
# ╟─4c85c2e9-5396-43d7-bea1-aed1414de316
