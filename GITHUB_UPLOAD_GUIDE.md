# Guía para publicar y ejecutar desde GitHub

## 1. Crear el repositorio

Cree un repositorio vacío en GitHub. No añada automáticamente README,
`.gitignore` ni licencia si va a cargar esta carpeta completa.

## 2. Subir desde RStudio

Abra `CalciumInsights.Rproj`, abra la terminal de RStudio y ejecute:

```bash
git init
git add .
git commit -m "Port current CalciumInsights app to golem"
git branch -M main
git remote add origin https://github.com/PROPIETARIO/REPOSITORIO.git
git push -u origin main
```

Sustituya `PROPIETARIO/REPOSITORIO` por el destino real.

## 3. Instalar desde otra laptop

En R o RStudio:

```r
install.packages("remotes")
remotes::install_github(
  "PROPIETARIO/REPOSITORIO",
  dependencies = TRUE,
  upgrade = "never"
)
CalciumInsights::run_app()
```

## 4. Trabajar con el código fuente en otra laptop

```bash
git clone https://github.com/PROPIETARIO/REPOSITORIO.git
cd REPOSITORIO
```

Abra `CalciumInsights.Rproj` y ejecute:

```r
source("install_dependencies.R")
shiny::runApp()
```

## 5. Comprobación antes de cada publicación

```r
devtools::document()
devtools::test()
devtools::check()
```

La carpeta `.github/workflows` contiene una comprobación automática para
pushes y pull requests. Una comprobación local exitosa debe realizarse antes
de depender del resultado de GitHub Actions.
