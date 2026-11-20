# Descripción del Modelo Reducido de Extinción de Apellidos

## Formulación General del Modelo

Consideremos un modelo dinámico de población a tiempo discreto $n \in \{0, 1, 2, \dots\}$, donde $n$ representa el número de generación. El modelo queda parametrizado por dos cantidades fundamentales:

$K \in \mathbb{N}$ ($K \ge 2$): Número inicial de apellidos distintos en la población, representados por el conjunto finito $\mathcal{A} = \{A_1, A_2, \dots, A_K\}$.

$N \in \mathbb{N}$**:** Tamaño total de la población, el cual se asume constante en cada generación (es decir, cada generación $n$ está conformada exactamente por $N$ individuos).

A diferencia de los modelos demográficos estándar (como el modelo multivariado de Wright-Fisher) que registran el vector completo de conteos e individuos por apellido $(N_1^{(n)}, \dots, N_K^{(n)} \text{ con } \sum N_i^{(n)} = N)$, el modelo reducido enfoca su atención únicamente en cuáles apellidos continúan vivos (es decir, aquellos que poseen al menos un representante en la población de la generación $n$).

Sea $X_n \subseteq \mathcal{A}$ el subconjunto de apellidos presentes en la generación $n$. El espacio de estados del proceso $(X_n)_{n \ge 0}$ está dado por el conjunto de partes no vacías de $\mathcal{A}$: $$\mathcal{E} = \{ S \subseteq \mathcal{A} : S \neq \emptyset \}$$

El tamaño del espacio de estados es exactamente $|\mathcal{E}| = 2^K - 1$. El conjunto vacío se excluye pues, por construcción, la población nunca se extingue totalmente ($N > 0$) y por ende al menos un apellido debe estar presente en cada generación.

## Supuesto de Herencia e Independencia (Hipótesis del Modelo Reducido)

Para que la secuencia del conjunto de apellidos vivos $(X_n)_{n \ge 0}$ constituya una Cadena de Markov homogénea por sí sola (sin necesidad de rastrear el número exacto de personas por apellido), se establece la siguiente regla de transmisión generacional:

Hipótesis Fundamental de Herencia Uniforme: Condicionado a que el conjunto de apellidos vivos en la generación $n$ sea $X_n = S$ con $|S| = k$ ($1 \le k \le K$), cada uno de los $N$ individuos de la generación $n+1$ elige su apellido de forma independiente y con distribución uniforme entre los $k$ apellidos disponibles en $S$.

Bajo este supuesto, la probabilidad de que un individuo cualquiera de la generación $n+1$ reciba un apellido específico $A_i \in S$ es exactamente: $$p = \frac{1}{|S|} = \frac{1}{k}$$

Debido a la independencia entre los $N$ individuos, la probabilidad condicional del nuevo conjunto de apellidos $X_{n+1}$ depende únicamente del conjunto de apellidos activos $X_n$ y no de las frecuencias relativas o proporciones de individuos que tenía cada apellido en la generación $n$, ni tampoco del historial de generaciones pasadas $X_0, X_1, \dots, X_{n-1}$. Esto satisface formalmente la propiedad de Markov:

$$\mathbb{P}(X_{n+1} = T \mid X_n = S, X_{n-1} = S_{n-1}, \dots, X_0 = S_0) = \mathbb{P}(X_{n+1} = T \mid X_n = S)$$

## Mecanismo de Inicialización (Generación 0)

La generación inicial $n = 0$ se construye bajo el mismo principio de distribución uniforme sobre los $K$ apellidos originales. Cada una de las $N$ personas de la generación $0$ recibe uno de los $K$ apellidos posibles de manera independiente con probabilidad $1/K$.

La probabilidad de que el estado inicial de la cadena sea un subconjunto específico $T \subseteq \mathcal{A}$ depende del número de personas $N$ y del número total de apellidos iniciales $K$.

# Probabilidades de Transición Combinatorias: La función $q(k, j)$

## Deducción Mediante el Principio de Inclusión-Exclusión

Sean $N$ personas que eligen su apellido de forma independiente y equiprobable entre $k$ apellidos disponibles. Queremos hallar la probabilidad de que el conjunto de apellidos que efectivamente aparecen en las $N$ personas sea exactamente un subconjunto fijado $T$ de tamaño $|T| = j$ (con $1 \le j \le k$).

El número total de asignaciones posibles de apellidos a la población es $k^N$. Para que el conjunto resultante sea exactamente $T$, debemos contar cuántas asignaciones asignan las $N$ personas únicamente a los $j$ apellidos de $T$ de tal forma que ningún apellido de $T$ quede sin al menos un representante (es decir, el número de funciones sobreyectivas de un conjunto de $N$ elementos en un conjunto de $j$ elementos).

Por el Principio de Inclusión-Exclusión, el número de funciones sobreyectivas $\operatorname{Sob}(N, j)$ está dado por:

$$\operatorname{Sob}(N, j) = \sum_{i=0}^{j} (-1)^{j-i} \binom{j}{i} i^N$$

Por lo tanto, la probabilidad de pasar a un estado específico $T$ de tamaño $j$ desde un estado de tamaño $k$ es:

$$q(k, j) = \frac{\operatorname{Sob}(N, j)}{k^N} = \frac{1}{k^N} \sum_{i=0}^{j} (-1)^{j-i} \binom{j}{i} i^N$$

Si estamos en un estado $S$ con $|S| = k$, la suma de las probabilidades hacia todos los subconjuntos posibles $T \subseteq S$ debe ser igual a 1:

$$\sum_{T \subseteq S} P(S, T) = \sum_{j=1}^{k} \binom{k}{j} q(k, j) = 1$$

## Distribución Inicial Generacional $\pi^{(0)}$

Asumiendo que la generación $n = 0$ se forma distribuyendo de manera uniforme las $N$ personas entre los $K$ apellidos originales, veamos que la probabilidad de que la cadena inicie en un estado específico $T \subseteq \mathcal{A}$ con $|T| = j$ viene dada directamente por evaluar $q(K, j)$:

$$\pi^{(0)}(T) = \mathbb{P}(X_0 = T) = q(K, |T|) = \frac{1}{K^N} \sum_{i=0}^{|T|} (-1)^{|T|-i} \binom{|T|}{i} i^N$$

En el instante inicial ($n = 0$), cada una de las $N$ personas de la población elige su apellido de manera independiente y equiprobable entre los $K$ apellidos disponibles. Puesto que cada uno de los $N$ individuos dispone de $K$ opciones independientes, el número total de asignaciones posibles de apellidos a la población es: 

$$\text{Casos totales} = K^N$$

Al ser las elecciones independientes y uniformes, las $K^N$ configuraciones poblacionales posibles son equiprobables, cada una con probabilidad de ocurrencia igual a $\frac{1}{K^N}$.

El estado de la cadena de Markov en el tiempo inicial, $X_0 = T$, representa el conjunto de apellidos vivos (aquellos que fueron elegidos por al menos una persona), donde $T \subseteq \mathcal{A}$ es un subconjunto no vacío de tamaño $|T| = j$ ($1 \le j \le K$).

Debido a que todos los $K$ apellidos iniciales tienen la misma probabilidad de ser seleccionados ($\frac{1}{K}$), el proceso posee una propiedad de simetría espacial:

La probabilidad $\pi^{(0)}(T) = \mathbb{P}(X_0 = T)$ no depende del contenido específico de los apellidos en $T$, sino únicamente de la cantidad de apellidos $j = |T|$ que contiene el subconjunto. Todos los subconjuntos del mismo tamaño $j$ poseen exactamente la misma probabilidad inicial.

Para determinar la probabilidad de que el conjunto de apellidos vivos en la generación inicial sea exactamente un subconjunto fijado $T$ de tamaño $|T| = j$, debemos contar cuántas de las $K^N$ asignaciones asignan a las $N$ personas exclusivamente a los $j$ apellidos de $T$, garantizando que ninguno de los $j$ apellidos quede sin al menos una persona. Esto equivale al número de funciones sobreyectivas de un conjunto de $N$ elementos (personas) a un conjunto de $j$ elementos (apellidos presentes). Por el Principio de Inclusión-Exclusión, este número viene dado por:

$$\operatorname{Sob}(N, j) = \sum_{i=0}^{j} (-1)^{j-i} \binom{j}{i} i^N$$

Al dividir entre el total de asignaciones posibles $K^N$, la probabilidad inicial de cualquier estado $T$ de tamaño $j$ es:

$$\pi^{(0)}(T) = \frac{\operatorname{Sob}(N, j)}{K^N} = \frac{1}{K^N} \sum_{i=0}^{j} (-1)^{j-i} \binom{j}{i} i^N$$

El comportamiento de $\pi^{(0)}$ depende fundamentalmente de la relación entre la población $N$ y el número de apellidos iniciales $K$. La probabilidad de que una persona dada no elija un apellido determinado es $\frac{K-1}{K}$. Por independencia, la probabilidad de que ninguna de las $N$ personas elija dicho apellido es: 

$$\mathbb{P}(\text{un apellido específico no aparezca}) = \left(\frac{K-1}{K}\right)^N = \left(1 - \frac{1}{K}\right)^N$$

Si la población $N$ es significativamente mayor que el número de apellidos $K$, el término $\left(1 - \frac{1}{K}\right)^N$ decrece exponencialmente hacia $0$. La probabilidad de que la generación inicial contenga el conjunto completo de todos los $K$ apellidos ($\pi^{(0)}(\mathcal{A})$) se aproxima a $1$: 

$$\pi^{(0)}(\mathcal{A}) \approx 1 - K \left(1 - \frac{1}{K}\right)^N$$

Las probabilidades de iniciar en estados con un número menor de apellidos ($j < K$) se vuelven extremadamente pequeñas.

Si el tamaño de la población es comparable o menor al número de apellidos, la probabilidad de que falten apellidos desde la generación inicial se vuelve considerablemente mayor, y la distribución inicial se dispersa sobre los subconjuntos de menor tamaño $j$.

# Construcción General de la Matriz de Transición $P$

Recordemos que el espacio de estados es $\mathcal{E} = \{ S \subseteq \mathcal{A} : S \neq \emptyset \}$, con $|\mathcal{E}| = 2^K - 1$. Dado que la población en la generación $n+1$ hereda sus apellidos únicamente de los apellidos que estaban vivos en la generación $n$ (estado $S$, con $|S| = k$), la probabilidad de transición $P(S, T) = \mathbb{P}(X_{n+1} = T \mid X_n = S)$ satisface: 

$$P(S, T) = \begin{cases} q(k, j) & \text{si } T \subseteq S, \\[4pt] 0 & \text{en otro caso} \end{cases}$$

donde $j = |T|$ representa la cantidad de apellidos vivos en el nuevo estado.

Si un apellido no está presente en $S$, ningún individuo de la generación $n+1$ puede heredarlo. Por lo tanto, no es posible que reaparezcan apellidos ni pasar a subconjuntos que no estén contenidos en $S$.

La probabilidad de que en la siguiente generación survivan exactamente los mismos $k$ apellidos es: 

$$P(S, S) = q(k, k) = \frac{\operatorname{Sob}(N, k)}{k^N}$$

Si en alguna generación la población queda reducida a un único apellido $S = \{A_i\}$ ($k = 1$), se cumple $q(1, 1) = 1$. Por ende, $P(\{A_i\}, \{A_i\}) = 1$, convirtiendo a los $K$ estados de tamaño 1 en estados absorbentes.

## Clasificación del Espacio de Estados

El espacio de estados $\mathcal{E}$ se particiona de manera natural en dos subconjuntos disjuntos:

Estados Transitorios ($\mathcal{T}$): Subconjuntos con 2 o más apellidos vivos ($|S| \ge 2$). Hay un total de $|\mathcal{T}| = 2^K - 1 - K$ estados transitorios. Desde cualquier estado transitorio existe una probabilidad positiva de perder apellidos a lo largo de las generaciones y nunca se puede regresar a un estado con mayor o igual número de apellidos.

Estados Absorbentes ($\mathcal{R}$): Subconjuntos de tamaño 1 ($|S| = 1$). Hay un total de $|\mathcal{R}| = K$ estados absorbentes (uno por cada apellido original). Una vez que la cadena entra en uno de estos estados, permanece en él indefinidamente.

## Forma Canónica

Para analizar el comportamiento asintótico de la cadena, ordenamos los estados de $\mathcal{E}$ colocando primero los $|\mathcal{T}|$ estados transitorios (en orden no creciente según su cardinalidad) y al final los $K$ estados absorbentes. Bajo esta permutación de filas y columnas, la matriz de transición de tamaño $(2^K - 1) \times (2^K - 1)$ adopta la Forma Canónica de Markov por Bloques:

$$P = \begin{pmatrix} Q & R \\[4pt] 0 & I_K \end{pmatrix}$$

donde:

$Q$ (Matriz de transiciones entre estados transitorios): De tamaño $(2^K - 1 - K) \times (2^K - 1 - K)$. Sus entradas $Q(S, T)$ representan la probabilidad de pasar del estado transitorio $S$ al estado transitorio $T$. Al ordenar los estados transitorios de mayor a menor cardinalidad, $Q$ resulta ser una matriz triangular superior por bloques. Sus elementos por debajo de la diagonal son estrictamente 0.

$R$ (Matriz de probabilidades de transición a la absorción): De tamaño $(2^K - 1 - K) \times K$. Sus entradas $R(S, \{A_i\})$ representan la probabilidad de pasar en un solo paso del estado transitorio $S$ al estado absorbente $\{A_i\}$.

$0$ (Matriz nula): De tamaño $K \times (2^K - 1 - K)$, refleja que es imposible salir de un estado absorbente hacia un estado transitorio.

$I_K$ (Matriz identidad): De tamaño $K \times K$, refleja que la probabilidad de permanecer en un estado absorbente es 1.

## La Matriz Fundamental $N$ (Tiempos Esperados en Estados Transitorios)

La matriz $N$, de tamaño $(2^K - 1 - K) \times (2^K - 1 - K)$, se define como la suma de la serie Neumann:

$$N = (I - Q)^{-1} = I + Q + Q^2 + Q^3 + \dots = \sum_{n=0}^{\infty} Q^n$$

La entrada $N_{S, T}$ representa el número promedio de generaciones que la población pasa en el estado transitorio $T$, habiendo comenzado en el estado $S$.

Sumando los elementos de la fila correspondiente al estado inicial $S$ en la matriz $N$, obtenemos el número esperado total de generaciones hasta que sobreviva un único apellido:

$$\tau(S) = \mathbb{E}[T_{\text{absorción}} \mid X_0 = S] = \sum_{T \in \mathcal{T}} N_{S, T}$$

## La Matriz de Probabilidades de Absorción $A$

La matriz $A$, de tamaño $(2^K - 1 - K) \times K$, se obtiene mediante el producto:

$$A = N R = (I - Q)^{-1} R$$

donde $A_{S, \{A_i\}}$` indica la probabilidad exacta de que la cadena termine en el estado absorbente $\{A_i\}$ (es decir, que el apellido $A_i$ sea el único superviviente) partiendo del estado transitorio $S$.

## Vector de Tiempos Esperados hasta la Absorción:

$$\tau = N \mathbf{1}$$

donde $\tau(S) = \sum_{T \in \mathcal{T}} N_{S, T}$ es la suma por filas de $N$, representando el tiempo total promedio que tarda la dinastía en quedar reducida a un único apellido.

### Cálculo Computacional

Para calcular estas matrices en un sistema computacional sin tener que invertir $(I - Q)$ directamente (lo cual puede ser inestable si hay problemas de condicionamiento), resolvemos los dos sistemas lineales equivalentes:

Para la Matriz Fundamental $N$:

$$(I - Q) N = I$$

Para la Matriz de Absorción $A$:

$$(I - Q) A = R$$

Dado que $Q$ es una matriz triangular superior al ordenar los estados por tamaño descendente (de 4 a 2 apellidos), la submatriz $I - Q$ también es triangular superior con unos en la diagonal principal. Esto significa que ambos sistemas se pueden resolver con sustitución hacia atrás (back-substitution) exacta, sin necesidad de eliminación gaussiana completa ni pivoteo.

# Clases de Comunicación e Irreducibilidad

Decimos que un estado $S$ alcanza a un estado $T$ ($S \to T$) si existe algún $n \ge 0$ tal que $P^n(S, T) > 0$. Decimos que $S$ y $T$ se comunican ($S \leftrightarrow T$) si $S \to T$ y $T \to S$.

Para dos estados transitorios $S, T \in \mathcal{T}$, si $T \subsetneq S$, la probabilidad de pasar en un solo paso es $P(S, T) = q(|S|, |T|) > 0$, por lo que $S \to T$. Una vez que un apellido se extingue, no puede regenerarse. Por ende, si $T \subsetneq S$, la probabilidad de regresar de $T$ a $S$ en cualquier número de pasos es estrictamente nula: $P^n(T, S) = 0$ para todo $n \ge 1$. Por lo tanto, $T \not\to S$. Si dos estados distintos tienen el mismo tamaño (por ejemplo, $S = ABC$ y $T = ABD$), es imposible pasar de $S$ a $T$ en ningún número de pasos ($P^n(S, T) = 0$), pues requeriría perder el apellido $C$ y ganar el apellido $D$ de forma simultánea.

Estados Transitorios ($2^K - 1 - K$ clases): Cada estado transitorio $S \in \mathcal{T}$ forma por sí solo una clase de comunicación no cerrada de un único elemento: $\{S\}$. No existen dos estados transitorios distintos que se comuniquen entre sí.

Estados Absorbentes ($K$ clases): Cada estado absorbente $\{A_i\} \in \mathcal{R}$ forma una clase de comunicación cerrada de un solo elemento: $\{\{A_i\}\}$, ya que $P(\{A_i\}, \{A_i\}) = 1$.

La cadena de Markov NO es irreducible (es una cadena reducible), ya que el espacio de estados $\mathcal{E}$ está fragmentado en $2^K - 1$ clases de comunicación independientes y no es posible alcanzar todos los estados desde cualquier punto de partida.

# Periodo de los Estados

El periodo de un estado $S \in \mathcal{E}$ se define como el máximo común divisor de los tiempos de retorno con probabilidad positiva:

$$d(S) = \gcd \{ n \ge 1 : P^n(S, S) > 0 \}$$

Estados Absorbentes: Para cada $\{A_i\} \in \mathcal{R}$, $P(\{A_i\}, \{A_i\}) = 1 > 0$, por lo que $d(\{A_i\}) = 1$.

Estados Transitorios: Para cualquier estado transitorio $S \in \mathcal{T}$ con $|S| = k \ge 2$, la probabilidad de permanecer en el mismo estado en un paso es:

$$P(S, S) = q(k, k) = \frac{\operatorname{Sob}(N, k)}{k^N} > 0$$

Al ser $P(S, S) > 0$, el número $n = 1$ pertenece al conjunto de tiempos de retorno, por lo cual $d(S) = 1$.

La cadena de Markov es aperiódica (todos sus estados tienen periodo $d = 1$).

# Matriz de Probabilidades de Llegada Eventual (Matriz $F$)

Definimos la entrada $F_{S, T}$ como la probabilidad de que la cadena visite eventualmente el estado $T$ al menos una vez, habiendo comenzado en el estado $S$: $$F_{S, T} = f_{S, T}^* = \mathbb{P}\left(\bigcup_{n=1}^{\infty} \{X_n = T\} \;\middle|\; X_0 = S\right)$$

$$F_{S, T} = f_{S, T}^* = \mathbb{P}\left(\bigcup_{n=1}^{\infty} \{X_n = T\} \;\middle|\; X_0 = S\right)$$

La matriz $F$ se relaciona directamente con la matriz fundamental $N = (I - Q)^{-1}$ y las probabilidades de absorción $A = NR$ mediante tres reglas teóricas:

### Transición entre estados transitorios distintos ($S \neq T \in \mathcal{T}$)

El número esperado de visitas $N_{S, T}$ equivale a la probabilidad de alcanzar $T$ por primera vez ($F_{S, T}$) multiplicada por el número promedio de visitas a $T$ partiendo de $T$ ($N_{T, T}$): 

$$N_{S, T} = F_{S, T} \cdot N_{T, T} \implies \mathbf{F_{S, T} = \frac{N_{S, T}}{N_{T, T}}}$$

### Retorno al mismo estado transitorio ($S = T \in \mathcal{T}$)

El número de visitas al propio estado inicial satisface la relación $N_{S, S} = \frac{1}{1 - f_{S, S}^*}$. Por lo tanto: 

$$f_{S, S}^* = 1 - \frac{1}{N_{S, S}} = P(S, S) = q(k, k)$$

Para cualquier estado $S \in \mathcal{T}$ (Transitorios), una vez que la cadena abandona el estado $S$ hacia un subconjunto propio $T \subsetneq S$, es imposible retornar a $S$.La probabilidad de retorno eventual a $S$ partiendo de $S$ es simplemente la probabilidad de no salir de él en el primer paso:

$$f_{S, S}^* = \mathbb{P}(\text{retornar a } S \mid X_0 = S) = P(S, S) = q(k, k) < 1$$

Como $f_{S, S}^* < 1$, por definición todos los estados con 2 o más apellidos son estrictamente transitorios.

### Llegada a estados absorbentes ($S \in \mathcal{T}, \{A_i\} \in \mathcal{R}$)

Como los estados absorbentes son trampas definitivas, la probabilidad de visitar el estado absorbente $\{A_i\}$ eventualmente coincide de forma exacta con la probabilidad de absorción $A_{S, \{A_i\}}$:

$$\mathbf{F_{S, \{A_i\}} = A_{S, \{A_i\}} = (NR)_{S, \{A_i\}}}$$

Un estado $S$ es recurrente si la probabilidad de regresar a él infinitas veces es 1. Para cada $\{A_i\} \in \mathcal{R}$, la probabilidad de retorno es $f_{\{A_i\}, \{A_i\}}^* = 1$. El tiempo medio de retorno es $\mathbb{E}[\tau_{\{A_i\}}^+] = 1 < \infty$, lo que los clasifica formalmente como recurrentes positivos.