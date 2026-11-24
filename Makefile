# Makefile - atajos para el proyecto Evolucion y Extincion de Apellidos.
#
# Multiplataforma via el script Julia subyacente. Requiere `julia` en el PATH.
# En Windows, si no tienes `make`, usa directamente los comandos `julia`
# documentados en el README (seccion "Notebooks de Pluto").
#
# Ejecuta `make` o `make help` para ver los objetivos disponibles.

JULIA ?= julia
PROJECT = --project=.

.DEFAULT_GOAL := help

.PHONY: help instantiate pluto bgw parejas bgw-dane parejas-dane experimentos simulation test test-markov bench

help: ## Muestra esta ayuda
	@echo "Objetivos disponibles:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

instantiate: ## Instala y precompila las dependencias del entorno (solo una vez)
	$(JULIA) $(PROJECT) -e "import Pkg; Pkg.instantiate()"

pluto: ## Lanza el servidor Pluto (elige el notebook en el navegador)
	$(JULIA) $(PROJECT) scripts/pluto.jl

bgw: ## Abre el notebook del modelo Bienayme-Galton-Watson
	$(JULIA) $(PROJECT) scripts/pluto.jl bgw

parejas: ## Abre el notebook basico del modelo por parejas (couple_model)
	$(JULIA) $(PROJECT) scripts/pluto.jl parejas

bgw-dane: ## Abre el notebook de resultados BGW con datos del Censo DANE 2018
	$(JULIA) $(PROJECT) scripts/pluto.jl bgw-dane

parejas-dane: ## Abre el notebook de resultados del arbol genealogico (DANE 2018)
	$(JULIA) $(PROJECT) scripts/pluto.jl parejas-dane

experimentos: ## Corre el catalogo de experimentos por parejas y guarda PNG + CSV en output/
	$(JULIA) $(PROJECT) scripts/experimentos_parejas.jl

simulation: ## Ejecuta la simulacion de demostracion (Galton-Watson)
	$(JULIA) $(PROJECT) scripts/run_simulation.jl

test: ## Corre los tests de los metodos numericos
	$(JULIA) $(PROJECT) test/runtests.jl

test-markov: ## Corre los tests del modelo de Markov / simulacion
	$(JULIA) $(PROJECT) test/test_markov.jl

bench: ## Ejecuta el benchmark de los metodos numericos
	$(JULIA) $(PROJECT) test/benchmarks.jl
