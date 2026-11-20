# Modelo Multiapellido reducido

## 1. Descripción del modelo

Consideramos una población con **4 apellidos** posibles, $\{A,B,C,D\}$, y **$N = 1000$ personas por generación**. Cada persona de la generación $n+1$ hereda su apellido de una persona de la generación $n$, así que la población se mantiene constante en 1000 personas.

Lo que nos interesa de cada generación **no es cuántas personas tiene cada apellido, sino cuáles apellidos siguen vivos** (es decir, aparecen al menos una vez).

### Supuesto de herencia (importante)

Para que el conjunto de apellidos vivos sea una cadena de Markov por sí sola, suponemos que:

> Cada una de las 1000 personas de la generación $n+1$ recibe su apellido **de forma independiente y uniforme entre los apellidos que siguen vivos en la generación $n$**.

Con este supuesto la probabilidad de transición depende solo de *qué apellidos están vivos* y no de cuántas personas tiene cada uno. Esto es lo que hace "reducido" al modelo. (Ver la sección 6 sobre qué pasa si se quita este supuesto.)

La generación 0 se genera igual, a partir de los 4 apellidos: cada persona recibe uno de los 4 con probabilidad $1/4$, de forma independiente.

## 2. Definición formal de la cadena

**Espacio de estados.** Sea $E$ la familia de subconjuntos no vacíos de $\{A,B,C,D\}$:

$$E = \{ S \subseteq \{A,B,C,D\} : S \neq \emptyset \}, \qquad |E| = 2^4 - 1 = 15.$$

El estado $S$ significa "en esa generación los apellidos presentes son exactamente los de $S$". Hay:

| Tamaño | Estados | Cantidad |
|---|---|---|
| 4 | ABCD | 1 |
| 3 | ABC, ABD, ACD, BCD | 4 |
| 2 | AB, AC, AD, BC, BD, CD | 6 |
| 1 | A, B, C, D | 4 |

El conjunto vacío no es un estado: siempre queda al menos un apellido.

**Proceso.** $X_n \in E$ es el conjunto de apellidos presentes en la generación $n$, $n = 0,1,2,\dots$

**Propiedades estructurales.**

- Un apellido extinto no reaparece: si $X_{n+1} = T$ y $X_n = S$, entonces $T \subseteq S$. Por eso, por ejemplo, $ABC \not\to ABCD$.
- Los estados de tamaño 1 ($A$, $B$, $C$, $D$) son **absorbentes**.
- Los 11 estados restantes (tamaño $\ge 2$) son **transitorios**: desde cada uno hay probabilidad positiva de perder apellidos y nunca se puede volver a un estado "más grande".

## 3. Probabilidad de que $m$ apellidos concretos cubran a las 1000 personas

Para calcular todo se necesita una única cantidad. Si 1000 personas reciben, de forma independiente y uniforme, un apellido entre $k$ disponibles, la probabilidad de que el conjunto de apellidos presentes sea **exactamente** un subconjunto concreto $T$ con $|T| = j$ es

$$q(k,j) = \frac{\operatorname{Sob}(N,j)}{k^{N}}, \qquad 
\operatorname{Sob}(N,j) = \sum_{i=0}^{j} (-1)^{j-i} \binom{j}{i}\, i^{N}.$$

$\operatorname{Sob}(N,j)$ es el número de funciones sobreyectivas de las $N$ personas sobre los $j$ apellidos de $T$ (por inclusión–exclusión): las asignaciones donde aparecen todos los apellidos de $T$ y ninguno fuera de $T$.

## 4. Distribución inicial $\pi^{(0)}$

Con 1000 personas y 4 apellidos equiprobables, para cada estado $T$:

$$\pi^{(0)}(T) = \mathbb{P}(X_0 = T) = q(4,|T|) = \frac{1}{4^{1000}} \sum_{i=0}^{|T|} (-1)^{|T|-i} \binom{|T|}{i}\, i^{1000}.$$

Por simetría, solo depende del tamaño de $T$:

| Estado(s) | Cuántos | $\pi^{(0)}$ de cada uno | Fórmula |
|---|---|---|---|
| ABCD | 1 | $1 - 4.605994\times 10^{-125}$ (aprox.) | $\dfrac{4^{1000} - 4\cdot 3^{1000} + 6\cdot 2^{1000} - 4}{4^{1000}}$ |
| ABC, ABD, ACD, BCD | 4 | $1.151499\times 10^{-125}$ | $\dfrac{3^{1000} - 3\cdot 2^{1000} + 3}{4^{1000}}$ |
| AB, AC, AD, BC, BD, CD | 6 | $9.332636\times 10^{-302}$ | $\dfrac{2^{1000} - 2}{4^{1000}}$ |
| A, B, C, D | 4 | $8.709810\times 10^{-603}$ | $\dfrac{1}{4^{1000}}$ |

Es decir, con 1000 personas es casi seguro que la generación inicial contiene los 4 apellidos. La probabilidad de que falte alguno es

$$1 - \pi^{(0)}(ABCD) = 4.605994\times 10^{-125} \approx 4\left(\tfrac{3}{4}\right)^{1000}.$$

Comprobación: la suma de las 15 probabilidades es exactamente 1.

## 5. Matriz de transición $P$

Si la generación $n$ tiene $k$ apellidos (estado $S$, $|S|=k$), la probabilidad de pasar a un estado $T$ es

$$P(S,T) = \begin{cases} q(k,|T|) & \text{si } T \subseteq S, \\ 0 & \text{en otro caso.} \end{cases}$$

Es decir, para pasar de $S$ a $T$ hace falta que las 1000 personas reciban apellidos de $T$ y que aparezcan todos los de $T$. En particular:

- **Quedarse igual** ($T = S$): todos los apellidos sobreviven, $P(S,S) = q(k,k)$.
- **Absorción**: $P(\{x\}, \{x\}) = q(1,1) = 1$.

Por ejemplo, $P(ABCD, ABC)$ es la probabilidad de que, habiendo los 4 apellidos en la generación $n$, en la $n+1$ no haya ninguno con apellido $D$ y aparezcan $A$, $B$ y $C$ al menos una vez: $q(4,3) = \dfrac{3^{1000} - 3\cdot 2^{1000} + 3}{4^{1000}}$.

### Valores de las probabilidades

Las entradas de $P$ solo dependen de $(k, j)$, el tamaño del estado de partida y el de llegada:

| Símbolo | Transición | Valor |
|---|---|---|
| $r_4$ | $ABCD \to ABCD$ | $1 - 4.605994\times 10^{-125}$ |
| $a$ | $ABCD \to$ cada triple | $1.151499\times 10^{-125}$ |
| $b$ | $ABCD \to$ cada par | $9.332636\times 10^{-302}$ |
| $c$ | $ABCD \to$ cada apellido solo | $8.709810\times 10^{-603}$ |
| $r_3$ | triple $\to$ el mismo triple | $1 - 2.431432\times 10^{-176}$ |
| $d$ | triple $\to$ cada par contenido | $8.104775\times 10^{-177}$ |
| $e$ | triple $\to$ cada apellido contenido | $7.563891\times 10^{-478}$ |
| $r_2$ | par $\to$ el mismo par | $1 - 1.866527\times 10^{-301}$ |
| $f$ | par $\to$ cada apellido contenido | $9.332636\times 10^{-302}$ |

(Cada fila suma 1; por ejemplo $r_4 + 4a + 6b + 4c = 1$ y $r_3 + 3d + 3e = 1$ exactamente.)

Los términos "$1-\varepsilon$" aparecen así porque $\varepsilon$ es tan pequeño que no cabe en un `double`: en coma flotante estándar $1 - 10^{-125}$ es exactamente $1$.

### Matriz completa (15 × 15)

Orden de los estados: ABCD, ABC, ABD, ACD, BCD, AB, AC, AD, BC, BD, CD, A, B, C, D.

| desde \ hacia | ABCD | ABC | ABD | ACD | BCD | AB | AC | AD | BC | BD | CD | A | B | C | D |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **ABCD** | $r_4$ | $a$ | $a$ | $a$ | $a$ | $b$ | $b$ | $b$ | $b$ | $b$ | $b$ | $c$ | $c$ | $c$ | $c$ |
| **ABC** | 0 | $r_3$ | 0 | 0 | 0 | $d$ | $d$ | 0 | $d$ | 0 | 0 | $e$ | $e$ | $e$ | 0 |
| **ABD** | 0 | 0 | $r_3$ | 0 | 0 | $d$ | 0 | $d$ | 0 | $d$ | 0 | $e$ | $e$ | 0 | $e$ |
| **ACD** | 0 | 0 | 0 | $r_3$ | 0 | 0 | $d$ | $d$ | 0 | 0 | $d$ | $e$ | 0 | $e$ | $e$ |
| **BCD** | 0 | 0 | 0 | 0 | $r_3$ | 0 | 0 | 0 | $d$ | $d$ | $d$ | 0 | $e$ | $e$ | $e$ |
| **AB** | 0 | 0 | 0 | 0 | 0 | $r_2$ | 0 | 0 | 0 | 0 | 0 | $f$ | $f$ | 0 | 0 |
| **AC** | 0 | 0 | 0 | 0 | 0 | 0 | $r_2$ | 0 | 0 | 0 | 0 | $f$ | 0 | $f$ | 0 |
| **AD** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $r_2$ | 0 | 0 | 0 | $f$ | 0 | 0 | $f$ |
| **BC** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $r_2$ | 0 | 0 | 0 | $f$ | $f$ | 0 |
| **BD** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $r_2$ | 0 | 0 | $f$ | 0 | $f$ |
| **CD** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $r_2$ | 0 | 0 | $f$ | $f$ |
| **A** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $1$ | 0 | 0 | 0 |
| **B** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $1$ | 0 | 0 |
| **C** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $1$ | 0 |
| **D** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | $1$ |

### Forma canónica

Ordenando los estados transitorios primero (los 11 de tamaño $\ge 2$) y los absorbentes después (A, B, C, D):

$$P = \begin{pmatrix} Q & R \\ 0 & I_4 \end{pmatrix}$$

con $Q$ de tamaño $11 \times 11$ (triangular, porque nunca se "gana" un apellido) y $R$ de tamaño $11 \times 4$. A partir de esto se calcula la matriz fundamental $M = (I - Q)^{-1}$ (tiempos esperados hasta la absorción) y las probabilidades de absorción $B = MR$ (qué apellido sobrevive).

## 6. Observaciones

1. **El supuesto de uniformidad es una simplificación.** En un modelo de herencia más realista (tipo Wright–Fisher), cada persona hereda el apellido de un padre elegido uniformemente *entre las personas*, y entonces la probabilidad de que un apellido desaparezca depende de **cuántas personas lo tienen**, no solo de si está presente. En ese caso el proceso "conjunto de apellidos vivos" deja de ser una cadena de Markov por sí solo (haría falta el vector de conteos). El modelo reducido es el que se obtiene si se supone que los apellidos vivos tienen la misma frecuencia en cada generación.
2. **Con $N=1000$ la extinción es casi imposible** en una sola generación (el término dominante es del orden de $10^{-125}$), así que la dinámica interesante aparece con $N$ pequeño o al revisar el modelo con frecuencias desiguales. Es útil dejar $N$ como parámetro.
3. **Observación numérica:** calcular $\operatorname{Sob}(N,j)/k^N$ con inclusión–exclusión directa en coma flotante sufre cancelación catastrófica y *overflow* ($4^{1000}\approx 10^{602}$). Los valores de este documento se calcularon con aritmética entera exacta; para $N$ grande en código conviene trabajar con enteros de precisión arbitraria o con logaritmos.
