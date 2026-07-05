# Notas de migración a golem

## Fuentes utilizadas

- Aplicación Shiny actual: `CalciumInsights_Modular_Fixed`.
- Estructura golem anterior: `CalciumInsights_AOG_lab`.

## Cambios de la aplicación

- La interfaz reproduce las dos pestañas activas de la app actual:
  `FFT + Baseline Analysis` y `Wavelet Ridgewalking`.
- La página Home usa el `index.html` y las imágenes de la app actual.
- La pestaña Home permanece oculta en la barra y la marca
  `CalciumInsights` funciona como botón de regreso.
- El módulo `Method Comparison` se conserva en el código, pero continúa
  desactivado, igual que en la app Shiny actual.
- Los módulos antiguos del paquete golem que no aparecen en la app actual
  no se copiaron al paquete actualizado.

## Adaptaciones necesarias para un paquete R

- Se eliminaron llamadas a `source()` dependientes del directorio de trabajo.
- Se eliminaron llamadas a `library()` ejecutadas al cargar archivos del
  paquete; las dependencias se declararon en `DESCRIPTION` y `NAMESPACE`.
- Se estableció `Collate` en `DESCRIPTION` para conservar el orden de carga:
  auxiliares, módulo FFT, módulo Wavelet, UI, servidor y lanzador.
- Los recursos web se movieron a `inst/app/www` y se exponen con
  `shiny::addResourcePath()`.
- El límite de carga de archivos de 200 MB se aplica al iniciar la app.
- Los datos experimentales del proyecto golem anterior se conservaron como
  ejemplos opcionales en `inst/extdata`.
- Se añadió un lanzador `app.R` para usar `shiny::runApp()` desde RStudio.
- Se añadió `run_app()` como función exportada para instalación desde GitHub.
- Se añadieron pruebas estructurales y un flujo de R CMD check para GitHub.

## Elementos no modificados

No se cambiaron fórmulas, umbrales, valores predeterminados, algoritmos de
detección, métricas, tablas descargables ni lógica de ajuste de los dos módulos
activos. Las únicas ediciones dentro de esos archivos fueron las necesarias
para sustituir `library()` y `source()` por la carga normal de un paquete.

## Elemento pendiente

La fuente golem anterior no contenía una licencia definitiva. El paquete usa
temporalmente `License: file LICENSE` con todos los derechos reservados. El
titular debe elegir una licencia antes de publicar el repositorio como proyecto
abierto.
