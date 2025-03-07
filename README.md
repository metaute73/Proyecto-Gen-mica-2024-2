# Proyecto-Gen-mica-2024-2
# Proyecto de Genómica para Clasificación de Pureza

## 📌 Descripción
Este proyecto desarrolla un software que clasifica experimentos de cultivo y secuenciación genómica en tres categorías:
- **Puro 🧫**: Una sola bacteria dominante sin contaminación significativa.
- **Mezcla de al menos dos 🧬🔬**: Presencia clara de dos o más bacterias.
- **Impuro 💀**: Una bacteria principal con restos significativos de otras secuencias.

## 🎯 Objetivo
Automatizar la clasificación de datos genómicos provenientes de experimentos de secuenciación, analizando métricas como longitud, profundidad, contenido GC y alelos alternativos.

## 📊 Datos de Entrada
El software procesa archivos **CSV** que contienen la siguiente información:
- **Scaffold**: Nombre del contig.
- **Length**: Longitud del contig.
- **GC**: Proporción de guanina y citosina.
- **Depth**: Profundidad media de lectura.
- **AltAllels**: Número de variantes detectadas en el contig.


## 🛠️ Funcionamiento
El software emplea métricas clave para evaluar la pureza del experimento:
1. **Métrica 1 - Dos Bacterias**: Penaliza diferencias significativas en longitudes y profundidades de los contigs.
2. **Métrica 2 - Contaminación**: Evalúa la cantidad de ADN ensamblado y la variabilidad genética.

Resultados Esperados
El software genera un reporte con la clasificación del experimento y métricas de respaldo.

Integrantes
Breadley Josshuen Marin Velasco
Daniel Metaute Medina
Daniela Hernandez Cardenas
Juan David Suarez Morales
Valentina Lopera Urrego

Cliente: Felipe Cabarcas Jaramillo

Este proyecto se distribuye bajo la licencia MIT.

Click here 👇🏼 to see the animations, the code and a detailed breakdown of everything!  
GD: <a href="https://colab.research.google.com/github/metaute73/Trabajo-01-RNA-corregido/blob/main/Punto_1_Optimizaci%C3%B3n_Num%C3%A9rica/Gradient_Descent___Evolution_Algorithms.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
