@echo off
title Pluto - Modelo Wright-Fisher
cd /d "%~dp0"

echo ========================================================
echo  Abriendo Notebook: Modelo Wright-Fisher
echo ========================================================
echo.

"%LOCALAPPDATA%\Programs\Julia-1.13.3\bin\julia.exe" --project=. -e "using Pluto; Pluto.run(notebook=joinpath(\"notebooks\", \"wright_fisher_model.jl\"))"

pause
