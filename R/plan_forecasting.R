##----------------------------------------------------------------
## Init by deleting all variables and functions
rm(list=ls())
## Set the working directory
setwd(".")
## To ensure reproducibility of the random numbers
set.seed(1)
## Packages used
## only needed in case you have not yet installed these packages
#update.packages(checkBuilt=TRUE, ask=FALSE)
library.path <- .libPaths()

print(library.path)
#install.packages(c("rugarch","rmgarch","xts","lubridate","LaplacesDemon","foreach", "doParallel"))#, dependencies = TRUE)
library(rugarch, quietly=TRUE, lib.loc = library.path)
library(rmgarch, quietly=TRUE,lib.loc = library.path)
library(xts, quietly=TRUE,lib.loc = library.path)
library(lubridate, quietly=TRUE, lib.loc = library.path)
library(LaplacesDemon, quietly=TRUE, lib.loc = library.path)
library(foreach, quietly=TRUE, lib.loc = library.path)
library(doParallel, quietly=TRUE, lib.loc = library.path)

## Source functions
sapply(dir("./R/plan functions",full.names=TRUE), source)

## Load the net-load data
df_load <- read.csv("./datasets/whole home power import 30 min 12 months averaged.csv", header=TRUE, row.names="timestamp")

## Drop NaNs
if (any(is.na(df_load))){
  df_load <- df_load[complete.cases(df_load), ]
}

## Load carbon intensity data
df_carbon <- read.csv("./datasets/Carbon_intensity_GB_2019.csv", header=TRUE, row.names="DATETIME")


## Load battery soc data
df_soc <- read.csv("./datasets/GB-EFR-30-m-2019-sonnenEco75.csv", header=TRUE,  row.names="dtm")
#rownames(df_soc) <- df$X

# convert to positive values
df_soc$soc_up = df_soc$soc_up*(-1)
df_soc$min_power = df_soc$min_power*(-1)

##----------------------------------------------------------------
##                            Forecasting
##----------------------------------------------------------------

T <- 48
refit <- 14
window_size <- 100
y_fit <- 21
r_type <- "I"
horizon <- 1
out_of_sample = (nrow(df_soc)/T) - window_size
out_of_sample

out_of_sample = (nrow(df_carbon)/T) - window_size
predict_quantiles(df_carbon, horizon, out_of_sample, T, window_size, refit, 100)
out_of_sample = (nrow(df_soc)/T) - window_size
predict_quantiles(df_soc, horizon, out_of_sample, T, window_size, refit, 100)
out_of_sample = (nrow(df_load)/T) - window_size
predict_quantiles(df_load, horizon, out_of_sample, T, window_size, refit, 100)
