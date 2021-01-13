##----------------------------------------------------------------
## Init by deleting all variables and functions
rm(list=ls())
## Set the working directory
setwd(".")
## To ensure reproducibility of the random numbers
set.seed(1) 
## Packages used
## only needed in case you have not yet installed these packages
#install.packages(c(""rugarch","rmgarch","xts","lubridate","LaplacesDemon","foreach", "doParallel"))  
library(rugarch)
library(rmgarch)
library(xts)
library(lubridate)
library(LaplacesDemon)
library(foreach)
library(doParallel)
#library(mvtnorm)
#library(VGAM)
#require(MASS)
#require(mvtnorm)
## Install GUROBI R-package 
#install.packages("doParallel") 


## Source functions
sapply(dir("plan functions",full.names=TRUE), source)


## Load the net-load data
df_load <- read.csv("whole home power import 30 min 12 months averaged.csv", header=TRUE, row.names="timestamp")
## Drop NaNs
if (any(is.na(df_load))){
  df_load <- df_load[complete.cases(df_load), ]
}
## Plot the example data 
plot(df_load$X25231, type="l")

df_load$X25231
dim(df_load)


## Load carbon intensity data
df_carbon <- read.csv("Carbon_intensity_GB_2019.csv", header=TRUE, row.names="DATETIME")#, nrows=8832)
tail(df_carbon)
dim(df_carbon)


## Load battery soc data
#df_soc <- read.csv("GB-EFR-30-m-2019-sonnenEco75.csv", header=TRUE,  row.names="dtm")#, nrows=8832) 
df_soc <- read.csv("GB-EFR-30-m-2019-tesla135.csv", header=TRUE,  row.names="dtm")#, nrows=8832) 
#rownames(df_soc) <- df$X
df_soc$soc_up = df_soc$soc_up*(-1)
df_soc$min_power = df_soc$min_power*(-1)
head(df_soc)
dim(df_soc)
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
predict_quantiles(df_load, horizon, out_of_sample, T, window_size, refit, 100) #[colnames(df_load)[1]]


# df_carbon_q <- predict_quantiles(df_carbon, horizon, out_of_sample, T, window_size, refit, 100)
# head(df_carbon_q)
# 
# df_soc_up_q <- predict_quantiles(df_soc[colnames(df_soc)[1]], horizon, out_of_sample, T, window_size, refit, 100)
# head(df_soc_up_q)
# 
# df_soc_down_q <- predict_quantiles(df_soc[colnames(df_soc)[2]], horizon, out_of_sample, T, window_size, refit, 100)
# head(df_soc_down_q)
# 
# df_pmax_q <- predict_quantiles(df_soc[colnames(df_soc)[3]], horizon, out_of_sample, T, window_size, refit, 100)
# head(df_pmax_q)
# 
# df_pmin_q <- predict_quantiles(df_soc[colnames(df_soc)[4]], horizon, out_of_sample, T, window_size, refit, 100)
# head(df_pmin_q)
# 
# df_load_q <- predict_quantiles(df_load[colnames(df_load)[2:4]], horizon, out_of_sample, T, window_size, refit, 100) #[colnames(df_load)[1]]
# head(df_load_q)


#matplot(df_load_q[49:(48+49),], type = "l",lwd=1, col = 'grey')


for (colname in colnames(df_load)[80]){
  print(colname)
  "Create dataframe for column"
  m <- matrix(df_load[[colname]], byrow = TRUE, nrow = nrow(df_soc)/T, ncol = T)
  rX = as.data.frame(m)
  rX['t'] <- seq(as.Date("2013-01-01 00:00:00"), by = "day", length.out = nrow(df_soc)/T)
  rX = xts(rX[, 1:T], order.by=as.Date(rX$t, format = "%Y-%m-%d", tz = "UTC"))
  "Predict using GARCH"
  dccrl <- predict_garch(rX, horizon, out_of_sample, T, window_size, refit)
  mu_matrix <- get_mu(dccrl, T)
  H_list <- get_cov(dccrl, T)
  m_real <- m[(window_size+1):dim(m)[1],]
  print(rmse(c(m_real)-c(mu_matrix)))
  
  "Create quantiles"
  quantiles <- seq(from=0.95, to=0.05, by=-0.05)
  m_quantiles <- matrix(, nrow = length(mu_matrix), ncol = length(quantiles)+1)
  m_quantiles[,length(quantiles)+1] <- matrix(t(mu_matrix), byrow=TRUE, nrow=length(mu_matrix), ncol=1)
  
  l_quantiles <- matrix(, nrow = length(mu_matrix), ncol = length(quantiles)+1)
  l_quantiles[,length(quantiles)+1] <- matrix(t(mu_matrix), byrow=TRUE, nrow=length(mu_matrix), ncol=1)
  for (i in 1:dim(m_real)[1]){
    U <- (chol((H_list[[i]][,,1])))
    mu <- matrix(c(rep.int(0,T)), nrow = T, ncol = 1, byrow = TRUE)
    X <- rmvnc(100, mu_matrix[i,], U)
    for (a in 1:length(quantiles)){
      m_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a], 
                                                      mean=colMeans(X), 
                                                      sd=c(apply(X, 2, sd)))
      
      l_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a],
                                                      mean=mu_matrix[i,],
                                                      sd=diag(chol(H_list[[i]][,,1])))
    }
  }
  "Create and save dataframe"
  df_q <- as.data.frame(m_quantiles)
  df_l <- as.data.frame(l_quantiles)
  colnames(df_q) <- append(quantiles, 'mean')
  #write.csv(df_q, sprintf('./forecasts/%s_forecast.csv',colname),row.names = FALSE)
}


cc <- c()
for (zz in 1:length(quantiles)){
  print(quantiles[zz])
  print(mean(pinLoss(matrix(m_real)[,1], df_q[,zz], quantiles[zz], add=FALSE)))
  cc <- append(cc, mean(pinLoss(matrix(m_real)[,1], df_q[,zz], quantiles[zz], add=FALSE)))
}
plot(cc)

cc <- c()
for (zz in 1:length(quantiles)){
  print(quantiles[zz])
  print(mean(pinLoss(matrix(m_real)[,1], df_l[,zz], quantiles[zz], add=FALSE)))
  cc <- append(cc, mean(pinLoss(matrix(m_real)[,1], df_l[,zz], quantiles[zz], add=FALSE)))
}
lines(cc)


for (i in 60){#1:dim(m_real)[1]){
  U <- (chol((H_list[[i]][,,1])))
  mu <- matrix(c(rep.int(0,T)), nrow = T, ncol = 1, byrow = TRUE)
  X <- rmvnc(100, mu_matrix[i,], U)
  matplot((t(X)), type = "l",lwd=1, col = 'grey')
  for (a in 1:length(quantiles)){
    m_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a], mean=mu_matrix[i,], sd=diag(chol(H_list[[i]][,,1])))
    lines(qnorm(p=quantiles[a], mean=mu_matrix[i,], sd=diag(chol(H_list[[i]][,,1]))), col='green', lwd=2)
    lines(qnorm(p=quantiles[a], mean=colMeans(X), sd=c(apply(X, 2, sd))), col='blue', lwd=2)
  }
  lines(m_real[i,], col='red',lwd=3)
  #lines(mu_matrix[i,], col='blue',lwd=3)
}


joint.density.plot(X[,1], X[,2], color=TRUE)

