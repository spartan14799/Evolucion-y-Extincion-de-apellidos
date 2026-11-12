@echo off
title Lanzador de Pluto - Evolucion y Extincion de Apellidos
cd /d "%~dp0"

echo ========================================================
echo  Iniciando Pluto.jl para los Notebooks del Proyecto...
echo  Se abrira una pestana en tu navegador automaticamente.
echo ========================================================
echo.

"%LOCALAPPDATA%\Programs\Julia-1.13.3\bin\julia.exe" --project=. -e "using Pluto; Pluto.run()"

pause
