##----------------------------------------------------------------
## Init by deleting all variables and functions
rm(list=ls())
## Set the working directory
setwd(".")
old_wd <- getwd()
## To ensure reproducibility of the random numbers
set.seed(1)

## Install GUROBI R-package
#install.packages("C:/gurobi902/win64/R/gurobi_9.0-2.zip", repos = NULL)

# only needed in case you have not yet installed these packages
install.packages(c("CVXR","matrixStats", "foreach", "doParallel"), quiet=TRUE)

## Packages used
library("CVXR", quietly=TRUE)
library("matrixStats", quietly=TRUE)
#library("gurobi", quietly=TRUE)
library("foreach", quietly=TRUE)
library("doParallel", quietly=TRUE)

## validate the installation
CVXR::installed_solvers()

## Source functions
sapply(dir("./R/plan functions",full.names=TRUE), source)

## Load the net-load data names
df_load <- read.csv("./datasets/load colnames.csv", header=TRUE, row.names="X")

bess <- "sonnen"
save_plans_dir <- file.path(getwd(),sprintf("plans_%s", bess))
save_sim_dir <- file.path(getwd(),sprintf("simulations_%s", bess))
#'%ni%' <- Negate('%in%')
soc_dir <- file.path(getwd(),sprintf("forecasts/%s", bess))

##--------------------------------------------------------------------------------------------#
##                                Battery optimization                                        #
##--------------------------------------------------------------------------------------------#

T <- 48
quantiles <- seq(from=0.95, to=0.05, by=-0.05)

## -------------------------------------- PRICES --------------------------------------------- ##

pb <- matrix(c(0.0652*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE) # 0.0285 0.057 - 10000 cycles/ 0.065 - sonnen
pfcr <- matrix(c(0.0978*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE)
pfit <- matrix(c(0.055*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE) # retailer tariff 5.5 pounds = 0.0524
ptou <- matrix(c(0.1662*rep.int(1,1), 0.1020*rep.int(1,14), 0.1662*rep.int(1,33)), nrow = T, ncol = 1, byrow = TRUE)

## --------------------------------- BATTERY STORAGE ---------------------------------------- ##

C <- matrix(c(3.3*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE) # 3.68 3.3
Q <- 0.9*7.5 # 13.5
eff <- 0.93
Q_min <- 0
Q_max <- Q
Q_start <- 0.5
Q_end <- 0.5

q <- Variable(T+1)
b = Variable(T)
bc = Variable(T)
bd = Variable(T)
b_ch_dis <- Variable(T, boolean=TRUE)

## ------------------------------------ NET LOAD -------------------------------------------- ##

nl <- Variable(T)
nl_neg = Variable(T)
nl_pos = Variable(T)
nl_pos_neg <- Variable(T, boolean=TRUE)
nl_fuse <- matrix(c(18.4*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE) # 80A or 18.4 kW, 13.2 - 60a

##----------------------------- ENHANED FREQUENCY RESPONSE ---------------------------------- ##

c_reg_share = Variable(T)
z = Variable(T)
z1 = Variable(T)

## --------------------------------- SUPPORTING ----------------------------------------- ##

e <- matrix(c(rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE)
dt <- 0.5

df_carbon_q <- read.csv(sprintf("./forecasts/%s_forecast.csv", 'CARBON_INTENSITY'), header=TRUE)
df_soc_up_q <- read.csv(file.path(soc_dir, sprintf("%s_forecast.csv", 'soc_up')), header=TRUE)
df_soc_down_q <- read.csv(file.path(soc_dir, sprintf("%s_forecast.csv", 'soc_down')), header=TRUE)
df_pmax_q <- read.csv(file.path(soc_dir, sprintf("%s_forecast.csv", 'max_power')), header=TRUE)
df_pmin_q <- read.csv(file.path(soc_dir, sprintf("%s_forecast.csv", 'min_power')), header=TRUE)

# uncomment to use multiple cores
#numCores <- detectCores()
#numCores
#registerDoParallel(numCores)  # use multicore, set to the number of our cores
#foreach (col_idx=1:dim(df_load)[2]) %dopar% {

foreach (col_idx=1:1){
  user <- colnames(df_load)[col_idx]
  df_load_q <- read.csv(sprintf("./forecasts/%s_forecast.csv", user), header=TRUE)
  for (day in 1:1){
  #for (day in 1:(dim(df_load_q)/T)[1]){
    print(c('User', col_idx, 'Day', day))
    dir.create(file.path(save_plans_dir, sprintf('%d',day)))
    user_results_iepos <- matrix(, nrow = (T+1), ncol = length(quantiles))
    dir.create(file.path(save_sim_dir, sprintf('%d',day)))
    user_results_online_sim <- matrix(, nrow = T, ncol = 5*length(quantiles))

    for (alpha in 1:length(quantiles)){
      ## ------------------------------------ NET LOAD ------------------------------------------- ##

      nl_alpha <- matrix(df_load_q[,alpha], nrow = nrow(df_load_q)/T,  ncol = T, byrow = TRUE)
      nl_alpha_day <- matrix(nl_alpha[day,], nrow = T,  ncol = 1, byrow = TRUE)
      nl_alpha_day[which(nl_alpha_day<(-4))] <- (-4)

      ## --------------------------------- CARBON INTENSITY ---------------------------------------- ##

      pco2_alpha <- matrix(df_carbon_q[,alpha], nrow = nrow(df_carbon_q)/T,  ncol = T, byrow = TRUE)
      pco2_alpha_day <- matrix(pco2_alpha[day,], nrow = T,  ncol = 1, byrow = TRUE)

      ##----------------------------- ENHANED FREQUENCY RESPONSE ------------------------------------##

      # down reg - neg power - positive energy - charge battery - withdraw power to the grid
      # up reg - pos power - neg energy - discharge battery - supply power to the grid

      max_p <- matrix(df_pmax_q[,alpha], nrow = nrow(df_pmax_q)/T,  ncol = T, byrow = TRUE)
      max_p[which(max_p<0)] <- 0
      C_reg_up <- matrix(max_p[day,], nrow = T,  ncol = 1, byrow = TRUE)

      min_p <- matrix(df_pmin_q[,alpha], nrow = nrow(df_pmin_q)/T,  ncol = T, byrow = TRUE)
      min_p[which(min_p<0)] <- 0
      C_reg_down <- matrix(min_p[day,], nrow = T,  ncol = 1, byrow = TRUE)

      Q_reg_delta_up <- matrix(df_soc_up_q[,alpha], nrow = nrow(df_soc_up_q)/T,  ncol = T, byrow = TRUE)
      Q_reg_delta_up_day <- matrix(Q_reg_delta_up[day,], nrow = T,  ncol = 1, byrow = TRUE)
      Q_reg_delta_up_day[which(Q_reg_delta_up_day<0)] <- 0

      Q_reg_delta_down <- matrix(df_soc_down_q[,alpha], nrow = nrow(df_soc_down_q)/T,  ncol = T, byrow = TRUE)
      Q_reg_delta_down_day <- matrix(Q_reg_delta_down[day,], nrow = T,  ncol = 1, byrow = TRUE)
      Q_reg_delta_down_day[which(Q_reg_delta_down_day<0)] <- 0

      ##------------------------------------ CONSTRAINS --------------------------------------------##

      constraints <- list(## Model the battery storage
                          q[1] == Q_start*Q,
                          q[T+1] == Q_end*Q, # $q_{t+1} = q_t + c_t$, $t = 1, \ldots, T − 1$,
                          diff(q) == -(bc*eff+bd/eff)*dt+c_reg_share*(Q_reg_delta_down_day-Q_reg_delta_up_day),
                          -q <= 0,#-Q_min*Q,
                          q <= Q,
                          bc <= 0,
                          -bd <= 0,
                          b==bc+bd,

                          ## Linearize binary and variable product for bd <= b_ch_dis*(C-c_reg_share*C_reg_up),
                          z >= 0,
                          z <= b_ch_dis,
                          -C_reg_up*(1-b_ch_dis) <= z - c_reg_share,
                          z - c_reg_share <= 0,
                          bd <= b_ch_dis*(C)-z*C_reg_up,

                          ## Linearize binary and variable product for bc >= (b_ch_dis-1)*(C-c_reg_share*C_reg_down),
                          z1 >= 0,
                          z1 <= b_ch_dis,
                          -(1-b_ch_dis) <= z1 - c_reg_share,
                          z1 - c_reg_share <= 0,
                          bc >= b_ch_dis*(C)-z1*C_reg_down-C+c_reg_share*C_reg_down,

                          ## Split net load in positive and negative for ToU and FiT charge
                          nl == nl_alpha_day - b + c_reg_share*(Q_reg_delta_down_day-Q_reg_delta_up_day)/dt,
                          nl == nl_pos + nl_neg,
                          nl_pos <= nl_pos_neg*nl_fuse, # limits on net-load power
                          -nl_neg <= -(nl_pos_neg-1)*nl_fuse, # limits on net-load power fuse size
                          -nl_pos <= 0,
                          nl_neg <= 0
      )

      const_fcr <- c(constraints, c(q[1:T] <= Q-c_reg_share*C*dt/2,
                                    q[1:T] >= c_reg_share*C*dt/2))

      ##------------------------------------ OBJECTIVES --------------------------------------------##

      obj_profit <- (## ToU
                      t(ptou)%*%(nl_pos-c_reg_share*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt +
                      ## Battery operation
                      t(pb)%*%(bd)*dt + t(pb)%*%(-bc)*dt + t(pb)%*%((Q_reg_delta_down_day+Q_reg_delta_up_day)*c_reg_share) +
                      ## FiT
                      t(pfit)%*%(nl_neg-c_reg_share*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt -
                      ## EFR
                      t(pfcr)%*%(c_reg_share*C)*dt +
                      t(ptou)%*%(Q_reg_delta_down_day*c_reg_share) + t(pfit)%*%(-Q_reg_delta_up_day*c_reg_share))

      obj_env <- t(pco2_alpha_day)%*%(nl)*dt
      obj_seff <- t(e)%*%(nl_pos-nl_neg)*dt #+ t(e)%*%(-q[1:T])
      objectives <- c(obj_profit, obj_env, obj_seff)

      dm_preferences <- c(0.273, 0.226, 0.501)
      n_obj <- length(objectives)
      weights <- diag(1, n_obj, n_obj)
      result_data <- matrix(, nrow = n_obj, ncol = n_obj)

      ##---------------------------- PARETO MINMAX APPROXIMATION -----------------------------------##
      for (i in 1:n_obj){
        if (i == 1){
          constr <- const_fcr
        }
        else{
          constr <- constraints
        }
        prob <- Problem(Minimize(weighted_sum(objectives, weights[i,])), constr)

        result <- solve(prob, solver = "SCS", verbose = FALSE, feastol = 1e-6, num_iter = Inf) #, feastol = 1e-6, num_iter = Inf)
        ### Fill in the objectives values
        result_data[i,1] <- ( ## ToU
                              t(ptou)%*%(result$getValue(nl_pos)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt +
                              ## Battery operation
                              t(pb)%*%(result$getValue(bd))*dt + t(pb)%*%(-result$getValue(bc))*dt +
                              t(pb)%*%((Q_reg_delta_down_day+Q_reg_delta_up_day)*result$getValue(c_reg_share)) +
                              ## FiT
                              t(pfit)%*%(result$getValue(nl_neg)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt -
                              ## EFR
                              t(pfcr)%*%(result$getValue(c_reg_share)*C)*dt +
                              t(ptou)%*%(Q_reg_delta_down_day*result$getValue(c_reg_share)) + t(pfit)%*%(-Q_reg_delta_up_day*result$getValue(c_reg_share)))

        result_data[i,2] <- t(pco2_alpha_day)%*%(result$getValue(nl))*dt
        result_data[i,3] <- t(e)%*%(result$getValue(nl_pos)-result$getValue(nl_neg))*dt  #+ t(e)%*%(-result$getValue(q)[1:T])
      }

      norm_objectives <- objectives
      for (obj_inx in 1:n_obj){
        norm_objectives[[obj_inx]] <- get_normalize_obj(result_data, objectives, dm_preferences, obj_inx)
      }

      prob <- Problem(Minimize(weighted_sum(norm_objectives, c(rep(1,n_obj)))), const_fcr)
      result <- solve(prob, solver = "SCS", verbose = FALSE, feastol = 1e-6, num_iter = Inf)

      gresult_data <- matrix(, nrow = 3, ncol = 1)
      gresult_data[1] <- (## ToU
                          t(ptou)%*%(result$getValue(nl_pos)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt +
                          ## Battery operation
                          t(pb)%*%(result$getValue(bd))*dt + t(pb)%*%(-result$getValue(bc))*dt +
                          t(pb)%*%((Q_reg_delta_down_day+Q_reg_delta_up_day)*result$getValue(c_reg_share)) +
                          ## FiT
                          t(pfit)%*%(result$getValue(nl_neg)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt -
                          ## EFR
                          t(pfcr)%*%(result$getValue(c_reg_share)*C)*dt +
                          t(ptou)%*%(Q_reg_delta_down_day*result$getValue(c_reg_share)) + t(pfit)%*%(-Q_reg_delta_up_day*result$getValue(c_reg_share)))

      gresult_data[2] <- t(pco2_alpha_day)%*%(result$getValue(nl))*dt
      gresult_data[3] <- t(e)%*%(result$getValue(nl_pos)-result$getValue(nl_neg))*dt #+ t(e)%*%(-result$getValue(q)[1:T])
      #gresult_data[4] <- t(e)%*%(-result$getValue(q)[1:T])

      user_results_iepos[1,alpha] <- result$value
      user_results_iepos[(2:(T+1)),alpha] <- result$getValue(nl)

      user_results_online_sim[(1:T),0*length(quantiles)+alpha] <- result$getValue(nl)
      user_results_online_sim[(1:T),1*length(quantiles)+alpha] <- result$getValue(b)
      user_results_online_sim[(1:T),2*length(quantiles)+alpha] <- result$getValue(c_reg_share)

      ## Flexibility availability

      free_reg_up <- C - result$getValue(c_reg_share)*C_reg_up - result$getValue(bd)
      energy_limit_up <- (eff*(result$getValue(q)[1:(T)] - result$getValue(c_reg_share)*C*0.5*dt))/dt
      flex_reg_up <- cbind(free_reg_up,energy_limit_up)
      flex_reg_up <- -matrixStats::rowMins(flex_reg_up)
      flex_reg_up

      free_reg_down <- -(C - result$getValue(c_reg_share)*C_reg_down - result$getValue(-bc))
      energy_limit_down <- (result$getValue(q)[1:(T)] - (Q-result$getValue(c_reg_share)*C*0.5*dt))/(eff*dt)
      flex_reg_down <- cbind(free_reg_down,energy_limit_down)
      flex_reg_down <- -matrixStats::rowMaxs(flex_reg_down)

      user_results_online_sim[(1:T),3*length(quantiles)+alpha] <- flex_reg_up
      user_results_online_sim[(1:T),4*length(quantiles)+alpha] <- flex_reg_down
    }

    setwd(sprintf('%s/%d/',save_plans_dir, day))
    for (l in 1:length(quantiles)){
     line<-paste(c(t(user_results_iepos)[l,1], paste(t(user_results_iepos)[l,2:(T+1)],collapse=",")), collapse = ':')
     write(line,file=sprintf('%s/%d/agent_%d.plans',save_plans_dir, day, grep(user, colnames(df_load))-1),append=TRUE)
    }

    user_timeseries <- as.data.frame(user_results_online_sim)
    colnames(user_timeseries) <- c(c(sprintf("nl_a%d", seq(1,length(quantiles)))),
                                   c(sprintf("pb_a%d", seq(1,length(quantiles)))),
                                   c(sprintf("c_reg_a%d", seq(1,length(quantiles)))),
                                   c(sprintf("flex_up_a%d", seq(1,length(quantiles)))),
                                   c(sprintf("flex_down_a%d", seq(1,length(quantiles)))))
    setwd(sprintf('%s/%d/',save_sim_dir, day))
    #print(getwd())
    write.csv(user_timeseries, sprintf('%s/%d/agent_%d.csv', save_sim_dir, day, grep(user, colnames(df_load))-1), row.names=FALSE)
    setwd(old_wd)
    #print(getwd())
  }
}
# uncomment to use multiple cores
# When you're done, clean up the cluster
#stopImplicitCluster()
