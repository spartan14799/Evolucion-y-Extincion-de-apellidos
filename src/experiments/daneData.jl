# src/experiments/daneData.jl
#
# DATOS reales del Censo DANE 2018 (Colombia): los 10 apellidos mas frecuentes,
# sus conteos y una paleta de colores. Archivo de solo datos, sin dependencia de
# ningun modelo; lo cargan los notebooks de notebooks/results/ y los runners
# daneBgwRunners.jl / daneCoupleRunners.jl.
#
# Documentacion completa (datos, uso): docs/code/daneData.md

using StatsBase

"""
    APELLIDOS_DANE :: Vector{String}

Los 10 apellidos mas frecuentes de Colombia (Censo Nacional de Poblacion y
Vivienda, DANE 2018), en orden de frecuencia descendente. Version con acentos.
"""
const APELLIDOS_DANE = [
    "Rodríguez", "Martínez", "García", "Gómez", "López",
    "Hernández", "Pérez", "Sánchez", "González", "Ramírez",
]

"""
    APELLIDOS_DANE_ASCII :: Vector{String}

Mismos apellidos que `APELLIDOS_DANE` sin acentos, para figuras cuya fuente no
renderiza bien los diacriticos (p. ej. el arbol genealogico).
"""
const APELLIDOS_DANE_ASCII = [
    "Rodriguez", "Martinez", "Garcia", "Gomez", "Lopez",
    "Hernandez", "Perez", "Sanchez", "Gonzalez", "Ramirez",
]

"""
    CONTEOS_DANE :: Vector{Int}

Numero de portadores de cada apellido de `APELLIDOS_DANE`, en el mismo orden,
segun el Censo DANE 2018.
"""
const CONTEOS_DANE = [
    754991, 602593, 574856, 569534, 565394,
    502626, 493764, 487112, 470505, 441179,
]

"""
    COLORES_DANE :: Vector{String}

Paleta de 10 colores (uno por apellido) para las figuras del arbol genealogico.
"""
const COLORES_DANE = [
    "#1E88E5", "#8E24AA", "#FB8C00", "#212121", "#6D4C41",
    "#E53935", "#43A047", "#00ACC1", "#FDD835", "#3949AB",
]

"""
    pesos_dane() -> Vector{Float64}

Peso relativo de cada apellido (conteo normalizado a suma 1), util para muestrear
poblaciones iniciales proporcionales al Censo DANE 2018.
"""
pesos_dane() = CONTEOS_DANE ./ sum(CONTEOS_DANE)

"""
    dane_dict() -> Dict{String,Int}

Diccionario apellido => conteo DANE 2018 (version con acentos).
"""
dane_dict() = Dict(zip(APELLIDOS_DANE, CONTEOS_DANE))
