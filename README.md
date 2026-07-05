# CalciumInsights — versión golem actualizada

CalciumInsights es una aplicación Shiny modular para el posprocesamiento,
visualización, detección de eventos y cuantificación de series temporales de
imágenes de calcio.

Esta versión convierte la app `CalciumInsights_Modular_Fixed` en un paquete
[golem](https://thinkr-open.github.io/golem/) listo para desarrollarse en
RStudio, almacenarse en GitHub e instalarse localmente como paquete de R.

## Módulos activos

1. **FFT + Baseline Analysis**
   - Suavizado FFT pasa-bajas.
   - Detección de picos.
   - Varias definiciones de línea base.
   - Análisis de sensibilidad a la línea base.
   - Métricas por evento y por traza.
   - AUC y modelos sigmoidales opcionales.
   - Tablas y gráficas descargables.

2. **Wavelet Ridgewalking**
   - Detección multiescala con wavelet Ricker/Mexican hat.
   - Construcción y filtrado de ridges.
   - Corrección de línea base posterior a la detección.
   - Métricas por evento y resumen.
   - AUC y modelos sigmoidales opcionales.
   - Tablas y gráficas descargables.

El archivo `R/mod_method_comparison.R` se conserva como en la app fuente, pero
el módulo permanece desactivado en la interfaz y el servidor.

## Requisitos

- R 4.1.0 o una versión posterior.
- RStudio Desktop es recomendado, aunque no obligatorio.
- Conexión a internet durante la primera instalación de dependencias.

## Ejecutar desde una carpeta descargada o clonada

1. Descargue o clone el repositorio.
2. Abra `CalciumInsights.Rproj` en RStudio.
3. En la consola de R ejecute:

```r
source("install_dependencies.R")
shiny::runApp()
```

`shiny::runApp()` utiliza el archivo `app.R` del repositorio, carga el paquete
en modo de desarrollo y abre la aplicación.

También puede ejecutar:

```r
source("dev/run_dev.R")
```

## Instalar directamente desde GitHub

Después de publicar el contenido de esta carpeta en un repositorio:

```r
install.packages("remotes")
remotes::install_github("AOG-Lab/CalciumInsights")
CalciumInsights::run_app()
```

Si usa otro propietario o nombre de repositorio, sustituya
`AOG-Lab/CalciumInsights` por `PROPIETARIO/REPOSITORIO`.

## Instalar el paquete desde una copia local

Desde la carpeta que contiene `DESCRIPTION`:

```r
install.packages("remotes")
remotes::install_local(".", dependencies = TRUE, upgrade = "never")
CalciumInsights::run_app()
```

## Formato de los datos

La aplicación acepta:

- `.csv`
- `.tsv`

La primera columna se interpreta como tiempo. Las columnas siguientes se
interpretan como señales de regiones de interés (ROI). La app no realiza
segmentación de imágenes ni extracción de ROI desde imágenes microscópicas.

Los archivos de ejemplo incluidos se encuentran en `inst/extdata`. Una vez
instalado el paquete, su ubicación puede consultarse con:

```r
system.file("extdata", package = "CalciumInsights")
```

## Estructura principal

```text
CalciumInsights/
├── app.R
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── app_ui.R
│   ├── app_server.R
│   ├── run_app.R
│   ├── mod_fft_baseline_sensitivity.R
│   ├── mod_wavelet_ridgewalking.R
│   └── utils_*.R
├── inst/
│   ├── app/www/
│   ├── extdata/
│   └── golem-config.yml
├── dev/run_dev.R
├── tests/testthat/
└── docs/
```

## Verificaciones recomendadas antes de publicar

Desde RStudio:

```r
devtools::document()
devtools::test()
devtools::check()
```

Luego pruebe ambos módulos con los datos simulados incluidos en la interfaz y
con al menos un archivo CSV real.

## Documentación adicional

Los documentos suministrados con la app actual se conservaron en `docs/`.
Consulte también:

- `MIGRATION_NOTES.md`
- `VALIDATION_REPORT.md`
- `FILE_MANIFEST.csv`

## Licencia

Los archivos suministrados no especificaban una licencia de código definitiva.
El archivo `LICENSE` conserva todos los derechos hasta que el titular seleccione
una licencia. Antes de distribuir públicamente el repositorio, sustituya ese
archivo y actualice el campo `License` de `DESCRIPTION`.
