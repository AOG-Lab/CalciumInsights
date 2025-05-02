# CalciumInsights

CalciumInsights is an interactive application built in R designed to analyze tissue-agnostic calcium traces.


![image](figures/CI_description.png)

# Installation

To install CalciumInsights R and RStudio are required:

<https://cran.r-project.org/>

<https://posit.co/download/rstudio-desktop/>

## Packages

1. **GOLEM**
2. **SHYNY**

## How to install the app

### Only the first time the app is installed, enter the following command in the console
```
install.packages("remotes")
```
### The following console command is to install from github
```
remotes::install_github("EMBRIOInstitute/CalciumInsights", auth_token = "your GitHub token")(FIXME)
```
## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CalciumInsights)
run_app()
```


