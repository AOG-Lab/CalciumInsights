# CalciumInsights

**CalciumInsights** is a modular R Shiny application for post-processing, visualization, event detection, baseline-sensitivity analysis, and quantification of calcium-imaging time-series data.

The application is distributed as a `golem` R package and can be installed directly on a Windows, macOS, or Linux laptop from the official GitHub repository:

**https://github.com/AOG-Lab/CalciumInsights**

## Active analysis modules

### FFT + Baseline Analysis

- FFT low-pass denoising
- Peak detection
- Multiple baseline definitions
- Baseline-sensitivity analysis
- Event-level and trace-level metrics
- Area-under-the-curve analysis
- Optional sigmoidal models
- Downloadable tables and plots

### Wavelet Ridgewalking

- Multiscale Ricker/Mexican-hat wavelet analysis
- Ridge construction and filtering
- Post-detection baseline correction
- Event-level metrics and summary statistics
- Area-under-the-curve analysis
- Optional sigmoidal models
- Downloadable tables and plots

`R/mod_method_comparison.R` is retained in the source code, but the method-comparison module is currently disabled in the user interface and server.

---

# 1. Install CalciumInsights directly from GitHub

This is the recommended installation method for most users. It installs CalciumInsights and its R package dependencies directly on the local laptop. The repository does not need to be downloaded or cloned manually.

## Step 1 — Install R

Install **R version 4.1.0 or later** from:

<https://cran.r-project.org/>

To check the installed R version, open R or RStudio and run:

```r
R.version.string
```

## Step 2 — Install RStudio Desktop

RStudio is recommended for running CalciumInsights locally:

<https://posit.co/download/rstudio-desktop/>

After installing R and RStudio, open **RStudio Desktop**.

## Step 3 — Install CalciumInsights

Copy the complete block below into the **RStudio Console** and press Enter:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never"
)
```

This command downloads the package source from the official GitHub repository and installs CalciumInsights and its dependencies into the local R library.

The installation may print many messages in the console. Warnings about packages being built are not necessarily errors. The installation is complete when the R prompt (`>`) appears again without an error message.

## Step 4 — Run CalciumInsights

After the installation is complete, run:

```r
CalciumInsights::run_app()
```

The application should open in the RStudio Viewer or in the default web browser.

To stop the application:

- Click **Stop** in RStudio, or
- Press the **Esc** key while the R console is active.

---

# 2. Open CalciumInsights after the first installation

The package only needs to be installed once. During later RStudio sessions, start the application with:

```r
CalciumInsights::run_app()
```

Alternatively:

```r
library(CalciumInsights)
run_app()
```

An internet connection is not normally required after the application and all dependencies have been installed, unless the application is updated or additional packages are needed.

---

# 3. Update CalciumInsights from GitHub

To replace an older local installation with the latest version available in the repository, run:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

Restart R after updating:

```text
RStudio > Session > Restart R
```

Then run:

```r
CalciumInsights::run_app()
```

To check the installed package version:

```r
packageVersion("CalciumInsights")
```

To locate the installed package on the laptop:

```r
find.package("CalciumInsights")
```

---

# 4. Remove CalciumInsights

To uninstall CalciumInsights from the local R library:

```r
remove.packages("CalciumInsights")
```

Restart RStudio after removing the package.

---

# 5. Required R packages

When CalciumInsights is installed with:

```r
remotes::install_github(
  "AOG-Lab/CalciumInsights",
  dependencies = TRUE
)
```

R installs the declared package dependencies automatically.

The application uses the following main R packages:

| Package | Purpose |
|---|---|
| `config` | Application configuration |
| `golem` | Shiny package framework |
| `shiny` | Web application framework |
| `shinyjs` | JavaScript utilities for Shiny |
| `DT` | Interactive data tables |
| `ggplot2` | Plotting |
| `dplyr` | Data manipulation |
| `vroom` | CSV and TSV import |
| `jsonlite` | JSON import |
| `pracma` | Peak detection and numerical functions |
| `prospectr` | Savitz–Golay filtering |
| `latex2exp` | Mathematical expressions in plots |
| `gridExtra` | Arrangement of multiple plots |
| `htmltools` | HTML components used by the interface |

The minimum versions explicitly required by the package include:

- R >= 4.1.0
- `config` >= 0.3.1
- `golem` >= 0.3.2
- `shiny` >= 1.7.2

Packages such as `pkgload`, `testthat`, and `devtools` are development tools. They are not required simply to start the installed application with `CalciumInsights::run_app()`.

---

# 6. Optional system build tools

Most Windows and macOS users receive precompiled R packages and do not need additional build tools. Install the appropriate tool only when R reports that a package must be compiled from source.

## Windows

Install the version of **Rtools** that matches the installed R version:

<https://cran.r-project.org/bin/windows/Rtools/>

Restart RStudio after installing Rtools.

To verify the installation:

```r
install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```

## macOS

Open the macOS Terminal and install the Apple Command Line Tools:

```bash
xcode-select --install
```

Restart RStudio after the installation finishes.

## Ubuntu or Debian Linux

Run the following commands in a terminal:

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

Make sure that the installed R version is 4.1.0 or later.

---

# 7. Alternative installation for developers

The direct GitHub installation described above is recommended for users who only need to run the application. Developers who want to inspect or modify the source code can clone or download the repository.

## Clone with Git

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

## Download the source as a ZIP file

1. Open <https://github.com/AOG-Lab/CalciumInsights>.
2. Select **Code**.
3. Select **Download ZIP**.
4. Extract the ZIP file.
5. Open `CalciumInsights.Rproj` in RStudio.
6. Run:

```r
source("install_dependencies.R")
shiny::runApp()
```

Do not run the project while it is still inside the compressed ZIP archive.

---

# 8. Input data format

The application accepts:

- `.csv`
- `.tsv`

The expected file structure is:

- The first column contains time values.
- Each additional column contains one calcium signal or region-of-interest trace.
- Column names should be unique.
- Analysis columns should contain numeric values.
- Missing or nonnumeric values should be reviewed before analysis.

Example:

```text
Time,ROI_1,ROI_2,ROI_3
0.00,0.102,0.118,0.095
0.10,0.106,0.121,0.099
0.20,0.111,0.125,0.101
```

CalciumInsights analyzes already extracted time-series traces. It does not segment microscopy images or extract regions of interest from raw microscopy images.

## Example data

After installing the package, locate the included example files with:

```r
system.file("extdata", package = "CalciumInsights")
```

List all example files with:

```r
list.files(
  system.file("extdata", package = "CalciumInsights"),
  full.names = TRUE
)
```

---

# 9. Troubleshooting

## Error: `there is no package called 'remotes'`

Install `remotes`:

```r
install.packages("remotes")
```

Then repeat:

```r
remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never"
)
```

## Error: `there is no package called 'pkgload'`

This error normally occurs when the application is being run from its source folder for development. Install `pkgload` with:

```r
install.packages("pkgload")
```

For a normal user installation, start the installed package instead:

```r
CalciumInsights::run_app()
```

## Error: `Failed to install 'CalciumInsights' from GitHub`

First restart RStudio and run:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("remotes")

remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

If the error states that build tools are missing, install Rtools on Windows or the Apple Command Line Tools on macOS as described above.

## Error: a dependency is missing

Install all declared dependencies again:

```r
remotes::install_github(
  repo = "AOG-Lab/CalciumInsights",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE
)
```

Or install a specific missing package:

```r
install.packages("PACKAGE_NAME", dependencies = TRUE)
```

## R cannot access CRAN

Set the CRAN repository explicitly:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
getOption("repos")
```

A firewall, proxy, institutional network, or antivirus configuration may block package downloads.

## GitHub rate-limit or authentication error

A public repository normally does not require a GitHub token. If GitHub reports a rate-limit error, create a GitHub personal access token and store it with the `gitcreds` package:

```r
install.packages("gitcreds")
gitcreds::gitcreds_set()
```

Then repeat the GitHub installation.

Do not place a personal access token directly in the README, source code, or a public repository.

## The browser does not open automatically

Run:

```r
CalciumInsights::run_app(launch.browser = TRUE)
```

When necessary, open the local address printed in the R console.

## Port already in use

Run the app on a different local port:

```r
CalciumInsights::run_app(port = 3839)
```

## The application opens but a data file is not accepted

Confirm that:

- The extension is `.csv` or `.tsv`.
- The first column contains time.
- The remaining analysis columns contain numeric values.
- Column names are unique.
- The delimiter matches the file content.

## Report installation information

When reporting an installation problem, include:

```r
sessionInfo()
packageVersion("CalciumInsights")
```

Also include the complete error message and the operating system being used.

---

# 10. Development checks

These steps are intended for developers and repository maintainers.

Install the development tools:

```r
install.packages(c("devtools", "testthat", "pkgload"))
```

Open `CalciumInsights.Rproj` and run:

```r
devtools::document()
devtools::test()
devtools::check()
```

Before publishing a new version:

1. Run the automated tests.
2. Run `devtools::check()`.
3. Test `CalciumInsights::run_app()` from a clean R session.
4. Test both active analysis modules.
5. Test at least one valid CSV or TSV file.
6. Confirm that tables and plots can be downloaded.

---

# 11. Repository

Official repository:

<https://github.com/AOG-Lab/CalciumInsights>

To report a problem or request a feature, use the repository's **Issues** section.

---

# 12. License

CalciumInsights is distributed under the MIT License. See the `LICENSE` and `LICENSE.md` files included in the repository for the applicable terms.
