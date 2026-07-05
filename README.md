# CalciumInsights

**CalciumInsights** is a modular R Shiny application for post-processing, visualization, event detection, baseline-sensitivity analysis, and quantification of calcium-imaging time-series data.

This repository is organized as a [`golem`](https://thinkr-open.github.io/golem/) R package. It can be opened directly in RStudio for development, run from the downloaded source folder, installed as a local R package, or installed from GitHub.

## Active analysis modules

### 1. FFT + Baseline Analysis

- FFT low-pass denoising
- Peak detection
- Multiple baseline definitions
- Baseline-sensitivity analysis
- Event-level and trace-level metrics
- Area-under-the-curve analysis
- Optional sigmoidal models
- Downloadable tables and plots

### 2. Wavelet Ridgewalking

- Multiscale Ricker/Mexican-hat wavelet analysis
- Ridge construction and filtering
- Post-detection baseline correction
- Event-level metrics and summary statistics
- Area-under-the-curve analysis
- Optional sigmoidal models
- Downloadable tables and plots

`R/mod_method_comparison.R` is retained from the source application, but the method-comparison module is currently disabled in the user interface and server.

---

# 1. System requirements

## Required software

Install the following software before installing the R packages:

1. **R version 4.1.0 or later**  
   Download R from the official CRAN website:  
   <https://cran.r-project.org/>

2. **RStudio Desktop** — recommended  
   Download the free desktop version from:  
   <https://posit.co/download/rstudio-desktop/>

3. **An internet connection**  
   It is needed during the first installation of R packages.

Git is optional. It is only required when cloning the repository with `git clone`. The project can also be downloaded from GitHub as a ZIP file.

No Python, Java, or MATLAB installation is required.

## Optional system build tools

Most Windows and macOS users receive precompiled R packages and will not need additional tools. Install the tools below only when R reports that a package must be compiled from source.

### Windows

Install the version of **Rtools** that matches your installed R version:

<https://cran.r-project.org/bin/windows/Rtools/>

After installing Rtools, restart RStudio and verify it with:

```r
install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```

### macOS

Install the Apple Command Line Tools from the macOS Terminal:

```bash
xcode-select --install
```

Restart RStudio after the installation finishes.

### Ubuntu or Debian Linux

Install R development tools and commonly required system libraries in a terminal:

```bash
sudo apt update
sudo apt install -y \
  r-base \
  r-base-dev \
  build-essential \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev
```

R 4.1.0 or later is required. On older Linux distributions, follow the current CRAN instructions for installing a recent R version.

---

# 2. R package dependencies

The project includes an installation script named `install_dependencies.R`. This is the recommended way to install all required packages.

## Packages required to run the application

| Package | Purpose |
|---|---|
| `config` | Application configuration |
| `golem` | Shiny package framework |
| `shiny` | Web application framework |
| `shinyjs` | JavaScript utilities for Shiny |
| `DT` | Interactive data tables |
| `ggplot2` | Plotting |
| `dplyr` | Data manipulation |
| `vroom` | Fast CSV and TSV import |
| `jsonlite` | JSON import |
| `pracma` | Peak detection and numerical functions |
| `prospectr` | Savitzky–Golay filtering |
| `latex2exp` | Mathematical expressions in plots |
| `gridExtra` | Arrangement of multiple plots |
| `htmltools` | HTML elements used by the interface |

The minimum versions explicitly required by the package are:

- R >= 4.1.0
- `config` >= 0.3.1
- `golem` >= 0.3.2
- `shiny` >= 1.7.2

## Packages used for local development, installation, and testing

| Package | Purpose |
|---|---|
| `pkgload` | Loads the package directly from its source folder |
| `testthat` | Runs automated tests |
| `remotes` | Installs the package locally or from GitHub |
| `httpuv` | Local Shiny web server support |

`devtools` is optional. Install it only when you plan to document, test, build, or check the package during development.

---

# 3. Recommended installation: download the project and run it in RStudio

This method is recommended for first-time users and for local development.

## Step 1 — Download and extract the project

From GitHub:

1. Open the repository page.
2. Select **Code**.
3. Select **Download ZIP**.
4. Extract the ZIP file to a normal folder on the laptop.

Do not run the project while it is still inside the compressed ZIP archive.

A suitable location is, for example:

```text
Windows: C:/Users/YourName/Documents/CalciumInsights
macOS:   /Users/YourName/Documents/CalciumInsights
Linux:   /home/YourName/CalciumInsights
```

Avoid placing the project in a temporary folder.

## Step 2 — Open the RStudio project

Open this file:

```text
CalciumInsights.Rproj
```

RStudio should open with the project folder as the working directory. In the **Files** pane, confirm that at least these files are visible:

```text
DESCRIPTION
app.R
install_dependencies.R
R/
inst/
```

You can verify the working directory in the R console:

```r
getwd()
file.exists("DESCRIPTION")
file.exists("app.R")
```

Both `file.exists()` calls should return `TRUE`.

## Step 3 — Select a CRAN repository

Run this command in the RStudio console:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
```

## Step 4 — Install all required R packages

Run:

```r
source("install_dependencies.R")
```

The script checks which packages are missing and installs only those packages.

It installs:

```r
c(
  "config",
  "golem",
  "shiny",
  "shinyjs",
  "DT",
  "ggplot2",
  "dplyr",
  "vroom",
  "jsonlite",
  "pracma",
  "prospectr",
  "latex2exp",
  "gridExtra",
  "htmltools",
  "pkgload",
  "testthat",
  "remotes",
  "httpuv"
)
```

Package installation is normally required only once for each R library.

## Step 5 — Run the application

Use either of the following methods.

### Method A: run through `app.R`

```r
shiny::runApp()
```

### Method B: run through the development launcher

```r
source("dev/run_dev.R")
```

The application should open in the RStudio Viewer or in the default web browser.

To stop the application, click **Stop** in RStudio or press the **Esc** key in the console.

---

# 4. Manual dependency installation

Use this method only when `source("install_dependencies.R")` does not complete successfully.

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(
  c(
    "config",
    "golem",
    "shiny",
    "shinyjs",
    "DT",
    "ggplot2",
    "dplyr",
    "vroom",
    "jsonlite",
    "pracma",
    "prospectr",
    "latex2exp",
    "gridExtra",
    "htmltools",
    "pkgload",
    "testthat",
    "remotes",
    "httpuv"
  ),
  dependencies = TRUE
)
```

Check whether any package is still missing:

```r
required_packages <- c(
  "config",
  "golem",
  "shiny",
  "shinyjs",
  "DT",
  "ggplot2",
  "dplyr",
  "vroom",
  "jsonlite",
  "pracma",
  "prospectr",
  "latex2exp",
  "gridExtra",
  "htmltools",
  "pkgload",
  "testthat",
  "remotes",
  "httpuv"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

missing_packages
```

An empty character vector means that all dependencies are installed:

```text
character(0)
```

---

# 5. Install CalciumInsights as a local R package

This method installs CalciumInsights into the user's R library. The app can then be started from any working directory.

Open `CalciumInsights.Rproj`, then run:

```r
install.packages("remotes")

remotes::install_local(
  path = ".",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

Start the installed application with:

```r
CalciumInsights::run_app()
```

To confirm that the package is installed:

```r
packageVersion("CalciumInsights")
find.package("CalciumInsights")
```

To remove the installed package:

```r
remove.packages("CalciumInsights")
```

---

# 6. Install CalciumInsights directly from GitHub

After the repository has been published on GitHub, install it with:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("remotes")

remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never"
)
```

Then run:

```r
CalciumInsights::run_app()
```

When the repository owner or repository name is different, replace:

```text
AOG-Lab/CalciumInsights
```

with:

```text
OWNER/REPOSITORY
```

## Update an existing GitHub installation

```r
remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

Restart R after updating the package, then run:

```r
CalciumInsights::run_app()
```

## Install a specific branch

```r
remotes::install_github(
  repo = "AOG-Lab/CalciumInsights@branch-name",
  dependencies = TRUE,
  upgrade = "never"
)
```

---

# 7. Clone the repository with Git

This method is optional.

Install Git from:

<https://git-scm.com/downloads>

Then run in a terminal:

```bash
git clone https://github.com/AOG-Lab/CalciumInsights.git
cd CalciumInsights
```

Open `CalciumInsights.Rproj` in RStudio and run:

```r
source("install_dependencies.R")
shiny::runApp()
```

---

# 8. Input data format

The application accepts:

- `.csv`
- `.tsv`

The expected structure is:

- The first column contains time values.
- Each additional column contains one calcium signal or region-of-interest trace.
- Column names should be unique.
- Values used for analysis should be numeric.
- Missing or nonnumeric values should be reviewed before analysis.

Example:

```text
Time,ROI_1,ROI_2,ROI_3
0.00,0.102,0.118,0.095
0.10,0.106,0.121,0.099
0.20,0.111,0.125,0.101
```

CalciumInsights analyzes already extracted time-series traces. It does not segment microscopy images or extract regions of interest from raw image files.

## Example data

Example files are included in:

```text
inst/extdata/
```

After installing the package, locate them with:

```r
system.file("extdata", package = "CalciumInsights")
```

List the example files with:

```r
list.files(
  system.file("extdata", package = "CalciumInsights"),
  full.names = TRUE
)
```

---

# 9. Project structure

```text
CalciumInsights/
├── app.R
├── CalciumInsights.Rproj
├── DESCRIPTION
├── NAMESPACE
├── install_dependencies.R
├── R/
│   ├── app_ui.R
│   ├── app_server.R
│   ├── run_app.R
│   ├── mod_fft_baseline_sensitivity.R
│   ├── mod_wavelet_ridgewalking.R
│   ├── mod_method_comparison.R
│   └── utils_*.R
├── inst/
│   ├── app/
│   │   └── www/
│   ├── extdata/
│   └── golem-config.yml
├── dev/
│   └── run_dev.R
├── tests/
│   └── testthat/
├── docs/
└── .github/
    └── workflows/
```

---

# 10. Development and package checks

These steps are optional and are intended for developers or repository maintainers.

Install `devtools`:

```r
install.packages("devtools")
```

Open `CalciumInsights.Rproj`, then run:

```r
devtools::document()
devtools::test()
devtools::check()
```

A standard local validation sequence is:

```r
pkgload::load_all()
testthat::test_dir("tests/testthat")
devtools::check()
```

Before publishing a new release:

1. Run the automated tests.
2. Run `devtools::check()`.
3. Start the app from the source folder.
4. Start the installed package with `CalciumInsights::run_app()`.
5. Test both active modules with included example data.
6. Test at least one valid CSV or TSV file from the intended workflow.
7. Confirm that tables and plots can be downloaded.

---

# 11. Troubleshooting

## Error: `there is no package called 'pkgload'`

Run:

```r
install.packages("pkgload")
```

Then:

```r
shiny::runApp()
```

## Error: one or more required packages are missing

Run:

```r
source("install_dependencies.R")
```

Or install the missing package directly:

```r
install.packages("PACKAGE_NAME", dependencies = TRUE)
```

## Error: `DESCRIPTION` or `app.R` cannot be found

The project was opened from the wrong folder.

Open `CalciumInsights.Rproj`, or set the working directory to the folder that contains `DESCRIPTION`:

```r
setwd("PATH/TO/CalciumInsights")
```

Then verify:

```r
file.exists("DESCRIPTION")
file.exists("app.R")
```

## Package installation fails on Windows

1. Update R to a supported version.
2. Restart RStudio.
3. Install the matching Rtools version if R reports that compilation tools are required.
4. Run:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
source("install_dependencies.R")
```

## Package installation fails on macOS

Install the Apple Command Line Tools:

```bash
xcode-select --install
```

Restart RStudio and repeat:

```r
source("install_dependencies.R")
```

## Package installation fails on Linux

Install the system requirements listed in the Linux section, restart R, and repeat:

```r
source("install_dependencies.R")
```

## R cannot access CRAN

Set the repository explicitly:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
```

Check the current setting:

```r
getOption("repos")
```

A firewall, proxy, institutional network, or antivirus configuration may also block downloads.

## The application opens but an uploaded file is not accepted

Confirm that:

- The file extension is `.csv` or `.tsv`.
- The first column contains time.
- The remaining analysis columns contain numeric values.
- Column names are unique.
- The delimiter matches the file extension and content.

## The browser does not open automatically

Run the source launcher, then open the local URL printed in the R console:

```r
shiny::runApp(launch.browser = TRUE)
```

## Port already in use

Run the app on another local port:

```r
shiny::runApp(port = 3839)
```

## Reset the R session

In RStudio, select:

```text
Session > Restart R
```

Then run:

```r
source("install_dependencies.R")
shiny::runApp()
```

## Report reproducibility information

Include the output of the following command when reporting an installation problem:

```r
sessionInfo()
```

Also include the complete error message and the operating system being used.

---

# 12. Additional documentation

Additional project documents are available in the repository:

- `MIGRATION_NOTES.md`
- `VALIDATION_REPORT.md`
- `GITHUB_UPLOAD_GUIDE.md`
- `FILE_MANIFEST.csv`
- `docs/`

---

# 13. License

The supplied source files did not specify a definitive open-source license. The current `LICENSE` file retains all rights until the copyright holder selects a distribution license.

Before publishing the repository publicly, update:

- `LICENSE`
- The `License` field in `DESCRIPTION`
- This README, when necessary
