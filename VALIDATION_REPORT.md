# Informe de validación del proyecto golem actualizado

## Resultado

Las comprobaciones estáticas realizadas sobre la carpeta final no detectaron
errores estructurales.

## Comprobaciones completadas

- Archivos R revisados: 24.
- Delimitadores y cadenas básicos: sin desequilibrios detectados.
- Campos obligatorios de `DESCRIPTION`: presentes.
- Campos duplicados en `DESCRIPTION`: ninguno.
- Archivos declarados en `Collate`: coinciden con todos los archivos de `R/`.
- Funciones principales de UI, servidor, FFT y Wavelet: presentes.
- Llamadas activas a `source()`, `setwd()`, `getwd()`, `library()` o
  `require()` dentro de `R/`: ninguna.
- Referencias `paquete::función`: cubiertas por `Imports` o por paquetes base.
- Recursos web copiados desde la app actual: hashes idénticos.
- Referencias locales de `index.html`: existentes.
- Archivos CSV de ejemplo legibles: 9 de 9.
- Metadatos de macOS (`__MACOSX`, `.DS_Store`, `._*`): excluidos del proyecto final.
- Repositorio `.git` de la app anterior: excluido.

## Integridad de la migración

Los dos archivos de módulos activos mantienen el contenido de la app Shiny
actual. Las diferencias se limitan a:

1. eliminar llamadas a `library()` del módulo FFT;
2. eliminar llamadas a `source()` dependientes del directorio de trabajo;
3. añadir comentarios que documentan la carga mediante `DESCRIPTION` y
   `NAMESPACE`.

La interfaz y el servidor golem llaman a:

- `mod_Denoising_data_fft_ui/server("fft_module")`;
- `mod_wavelet_ridgewalking_ui/server("wavelet_module")`.

El módulo de comparación continúa inactivo, igual que en la app fuente.

## Datos de ejemplo

| Archivo | Filas | Columnas |
|---|---:|---:|
| `Arabidopsis1.csv` | 361 | 6 |
| `Arabidopsis2.csv` | 361 | 6 |
| `Arabidopsis3.csv` | 361 | 6 |
| `Drosophila1.csv` | 361 | 11 |
| `Drosophila2.csv` | 361 | 11 |
| `Drosophila3.csv` | 361 | 11 |
| `Zebrafish1.csv` | 90 | 6 |
| `Zebrafish2.csv` | 90 | 6 |
| `Zebrafish3.csv` | 90 | 6 |

## Limitación de esta validación

El entorno utilizado para preparar el proyecto no contiene una instalación
ejecutable de R. Por esa razón no fue posible iniciar una sesión Shiny ni
ejecutar `R CMD check`, `devtools::test()` o pruebas interactivas de los
controles y descargas.

Después de abrir el proyecto en una computadora con R, ejecute:

```r
source("install_dependencies.R")
devtools::document()
devtools::test()
devtools::check()
shiny::runApp()
```

La carpeta `.github/workflows` también ejecutará R CMD check cuando el proyecto
se publique en GitHub.
