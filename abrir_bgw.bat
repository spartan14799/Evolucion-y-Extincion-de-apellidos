@echo off
title Pluto - Modelo Bienayme-Galton-Watson
cd /d "%~dp0"

echo ========================================================
echo  Abriendo Notebook: Modelo Bienayme-Galton-Watson
echo ========================================================
echo.

"%LOCALAPPDATA%\Programs\Julia-1.13.3\bin\julia.exe" --project=. -e "using Pluto; Pluto.run(notebook=joinpath(\"notebooks\", \"bgw_model.jl\"))"

pause
