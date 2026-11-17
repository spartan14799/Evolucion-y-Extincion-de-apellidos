# scripts/pluto.jl
#
# Lanzador multiplataforma de Pluto.jl para los notebooks del proyecto.
# Reemplaza a los antiguos .bat (que ataban la ruta de Julia a Windows).
#
# Uso (desde la raiz del proyecto):
#     julia --project=. scripts/pluto.jl          # servidor Pluto (elige notebook en el navegador)
#     julia --project=. scripts/pluto.jl bgw      # abre notebooks/bgw_model.jl
#     julia --project=. scripts/pluto.jl wf        # abre notebooks/wright_fisher_model.jl
#
# No hardcodea la ruta de Julia: usa el interprete con el que se lanza y el
# entorno activado por --project=. Antes de abrir instancia el entorno si hace
# falta, para que funcione en una maquina recien clonada.

import Pkg

# Instancia el entorno del proyecto si todavia no esta resuelto/precompilado.
# Pkg.instantiate() es idempotente: si ya esta listo, no hace trabajo extra.
Pkg.instantiate()

using Pluto

const NOTEBOOKS = Dict(
    "bgw" => joinpath("notebooks", "bgw_model.jl"),
    "wf"  => joinpath("notebooks", "wright_fisher_model.jl"),
)

function main(args)
    if isempty(args)
        # Sin argumento: servidor Pluto general, el usuario elige en el navegador.
        Pluto.run()
        return
    end

    clave = lowercase(args[1])
    if !haskey(NOTEBOOKS, clave)
        validas = join(sort(collect(keys(NOTEBOOKS))), ", ")
        println(stderr, "Notebook desconocido: '$(args[1])'. Opciones validas: $validas")
        exit(1)
    end

    ruta = NOTEBOOKS[clave]
    if !isfile(ruta)
        println(stderr, "No se encontro el notebook en la ruta esperada: $ruta")
        println(stderr, "Ejecuta el comando desde la raiz del proyecto.")
        exit(1)
    end

    Pluto.run(notebook = ruta)
end

main(ARGS)
