# Prior to this installation, make sure your R version is up to date!! 

# Install INLA ####

install.packages("INLA", 
                 repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), 
                 dep=TRUE)


# Install inlabru ####

# Enable universe(s) by inlabru-org
options(repos = c(
  inlabruorg = "https://inlabru-org.r-universe.dev",
  INLA = "https://inla.r-inla-download.org/R/testing",
  CRAN = "https://cloud.r-project.org"
))

# Install package
install.packages("inlabru")

# Install additional needed packages ####

install.packages(c(
  "CARBayesdata",
  "dplyr",
  "fmesher",
  "ggplot2",
  "lubridate",
  "mapview",
  "patchwork",
  "scico",
  "sdmTMB",
  "sf",
  "spatstat",
  "spdep",
  "terra",
  "tidyr",
  "tidyterra",
  "viridis",
  "palmerpenguins"
))


# Installation check ####

## INLA ####
df <- data.frame(y = rnorm(100) + 10)

fit <- INLA::inla(
  y ~ 1,
  data = df
)

summary(fit)

## inlabru ####

fit <- inlabru::bru(
  y ~ Intercept(1, prec.linear = exp(-7)),
  data = df
)

summary(fit)
