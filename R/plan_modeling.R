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
#install.packages(c("CVXR","matrixStats", "foreach", "doParallel"), quiet=TRUE)   # only needed in case you have not yet installed these packages

## Packages used
library(CVXR, quietly=TRUE)
library(matrixStats, quietly=TRUE)
library(gurobi, quietly=TRUE)
library(foreach, quietly=TRUE)
library(doParallel, quietly=TRUE)

## validate the installation
CVXR::installed_solvers()

## Source functions
sapply(dir("plan functions",full.names=TRUE), source)

## Load the net-load data
df_load <- read.csv("load colnames.csv", header=TRUE, row.names="X")
## Drop NaNs
#if (any(is.na(df_load))){
#  df_load <- df_load[complete.cases(df_load), ]
#}


bess <- "sonnen"
save_plans_dir <- file.path(getwd(),sprintf("plans_%s", bess))
save_sim_dir <- file.path(getwd(),sprintf("simulations_%s", bess))
#'%ni%' <- Negate('%in%')
soc_dir <- file.path(getwd(),sprintf("forecasts/%s", bess))

#print(soc_dir)
#file.access(".", 2)

##--------------------------------------------------------------------------------------------#
##                                Battery optimization                                        #
##-------------------------------------------------------------------------------------------- #

T <- 48
quantiles <- seq(from=0.95, to=0.05, by=-0.05)

## -------------------------------------- PRICES --------------------------------------------- ##

pb <- matrix(c(0.0652*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE) #0.0285 0.057 - 10000 cycles/ 0.065 - sonnen
pfcr <- matrix(c(0.0978*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE)
pfit <- matrix(c(0.055*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE) # octopus 5.5 # 0.0524
ptou <- matrix(c(0.1662*rep.int(1,1), 0.1020*rep.int(1,14), 0.1662*rep.int(1,33)), nrow = T, ncol = 1, byrow = TRUE) #16.62 10.20 octopus / Npower 20.19 , 10.61

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

##----------------------------- ENHANED FREQUENCY RESPONSE ------------------------------------##

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

#print(head(df_pmin_q))
#print(file.path(soc_dir, sprintf("%s_forecast.csv", 'soc_up')))
#stopImplicitCluster()

#system.time({
#X95046 - 4kW
numCores <- detectCores()
numCores
registerDoParallel(numCores)  # use multicore, set to the number of our cores

foreach (col_idx=1:dim(df_load)[2]) %dopar% {  #1:dim(df_load)[2] , .export=ls(envir=globalenv()) #:dim(df_load)[2]
#for (col_idx in 1){                                                           #dim(df_load)[2]
  user <- colnames(df_load)[col_idx]
  #print(user)
  #print(col_idx)
  df_load_q <- read.csv(sprintf("./forecasts/%s_forecast.csv", user), header=TRUE)
  #print(head(df_load_q))
  for (day in 101:150){ #(dim(df_load_q)/T)[1]){
  #foreach(day=2:10, .export=ls(envir=globalenv())) %dopar% {                                      #1:(dim(df_load_q)/T)[1]){
    print(c('User', col_idx, 'Day', day))
    #print(file.path(save_plans_dir, sprintf('%d',day)))
    dir.create(file.path(save_plans_dir, sprintf('%d',day)))
    user_results_iepos <- matrix(, nrow = (T+1), ncol = length(quantiles))
    #print(user_results_iepos)
    dir.create(file.path(save_sim_dir, sprintf('%d',day)))
    user_results_online_sim <- matrix(, nrow = T, ncol = 5*length(quantiles))

    for (alpha in 1:length(quantiles)){
      #print(alpha)
      ## ------------------------------------ NET LOAD ------------------------------------------- ##

      nl_alpha <- matrix(df_load_q[,alpha], nrow = nrow(df_load_q)/T,  ncol = T, byrow = TRUE)
      nl_alpha_day <- matrix(nl_alpha[day,], nrow = T,  ncol = 1, byrow = TRUE)
      nl_alpha_day[which(nl_alpha_day<(-4))] <- (-4)

      #load_alpha_day <- nl_alpha_day
      #solar_alpha_day <- nl_alpha_day

      #load_alpha_day[which(load_alpha_day<0)] <- 0
      #solar_alpha_day[which(solar_alpha_day>0)] <- 0

      ## --------------------------------- CARBON INTENSITY ---------------------------------------- ##

      pco2_alpha <- matrix(df_carbon_q[,alpha], nrow = nrow(df_carbon_q)/T,  ncol = T, byrow = TRUE)
      pco2_alpha_day <- matrix(pco2_alpha[day,], nrow = T,  ncol = 1, byrow = TRUE)

      ##----------------------------- ENHANED FREQUENCY RESPONSE ------------------------------------##

      #down reg - neg power - positive energy - charge battery - withdraw power to the grid
      #up reg - pos power - neg energy - discharge battery - supply power to the grid

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
      #obj_blackout <- t(e)%*%(-q[1:T])
      objectives <- c(obj_profit, obj_env, obj_seff)#, obj_blackout)

      #dm_preferences <- c(0.468+0.01625, 0.216+0.01625, 0.236+0.01625, 0.015+0.01625)
      dm_preferences <- c(0.273, 0.226, 0.501)
      #dm_preferences <- c(0.99, 0., 0., 0.)
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

        result <- solve(prob, solver = "GUROBI", verbose = FALSE, feastol = 1e-6, num_iter = Inf) #, feastol = 1e-6, num_iter = Inf)
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
        #result_data[i,4] <- t(e)%*%(-result$getValue(q)[1:T])
      }

      norm_objectives <- objectives
      for (obj_inx in 1:n_obj){
        norm_objectives[[obj_inx]] <- get_normalize_obj(result_data, objectives, dm_preferences, obj_inx)
      }

      prob <- Problem(Minimize(weighted_sum(norm_objectives, c(rep(1,n_obj)))), const_fcr) #constraints) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check !!!!!!!!!!!!!!!!!!!!!!!
      result <- solve(prob, solver = "GUROBI", verbose = FALSE, feastol = 1e-6, num_iter = Inf)
      #if (result$status %ni% "optimal"){
      #  print('Status is not optimal')
      #  print(c('User', col_idx, 'Day', day))
      #}
      #result$status
      #result$solve_time
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
    #print(getwd())
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
# When you're done, clean up the cluster
stopImplicitCluster()


#user_timeseries <- read.csv(sprintf('./%s/%d/agent_%d.csv', "simulations", 15, 10))
#user_results_iepos <- table(sprintf('./%s/%d/agent_%d.plans', "plans", 15, 100), sep = ":", header = NULL)
#user_results_iepos
#
#matplot((user_results_iepos)[2:(T+1),], type="l")
#matplot((user_results_iepos)[1,], type="l")
#
#
#user_timeseries
#matplot(user_timeseries[,1:19], type="l")
#matplot(user_timeseries[,20:38], type="l")
#matplot(user_timeseries[,39:57], type="l")
#matplot(user_timeseries[,58:76], type="l")
#matplot(user_timeseries[,87:95], type="l")
#plot(result$getValue(q)[1:(T)])


# i <- 10
# alpha <- 10 ## 1 (0.05) - the most risk, 10 (0.95) - the most conservative


# ## ------------------------------------ NET LOAD ------------------------------------------- ##
#
# # df_load_q <- read.csv("./forecasts/X95141_forecast.csv", header=TRUE)
# # head(df_load_q)
#
# nl_alpha <- matrix(df_load_q[,alpha], nrow = nrow(df_load_q)/T,  ncol = T, byrow = TRUE)
# nl_alpha_day <- matrix(nl_alpha[i,], nrow = T,  ncol = 1, byrow = TRUE)
# load_alpha_day <- nl_alpha_day
# solar_alpha_day <- nl_alpha_day
#
# load_alpha_day[which(load_alpha_day<0)] <- 0
# solar_alpha_day[which(solar_alpha_day>0)] <- 0
# solar_alpha_day[which(solar_alpha_day<(-2))] <- (-2)
#
#
#
# ## --------------------------------- CARBON INTENSITY ---------------------------------------- ##
#
# pco2_alpha <- matrix(df_carbon_q[,alpha], nrow = nrow(df_carbon_q)/T,  ncol = T, byrow = TRUE)
# pco2_alpha_day <- matrix(pco2_alpha[i,], nrow = T,  ncol = 1, byrow = TRUE)
# pco2_alpha_day
#
#
# ##----------------------------- ENHANED FREQUENCY RESPONSE ------------------------------------##
#
# #down reg - neg power - positive energy - charge battery - withdraw power
# #up reg - pos power - neg energy - discharge battery - supply power
#
# max_p <- matrix(df_pmax_q[,alpha], nrow = nrow(df_pmax_q)/T,  ncol = T, byrow = TRUE)
# max_p[which(max_p<0)] <- 0
# C_reg_up <- matrix(max_p[i,], nrow = T,  ncol = 1, byrow = TRUE)
#
# min_p <- matrix(df_pmin_q[,alpha], nrow = nrow(df_pmin_q)/T,  ncol = T, byrow = TRUE)
# min_p[which(min_p<0)] <- 0
# C_reg_down <- matrix(min_p[i,], nrow = T,  ncol = 1, byrow = TRUE)
#
# Q_reg_delta_up <- matrix(df_soc_up_q[,alpha], nrow = nrow(df_soc_up_q)/T,  ncol = T, byrow = TRUE)
# Q_reg_delta_up_day <- matrix(Q_reg_delta_up[i,], nrow = T,  ncol = 1, byrow = TRUE)
# Q_reg_delta_up_day[which(Q_reg_delta_up_day<0)] <- 0
#
# Q_reg_delta_down <- matrix(df_soc_down_q[,alpha], nrow = nrow(df_soc_down_q)/T,  ncol = T, byrow = TRUE)
# Q_reg_delta_down_day <- matrix(Q_reg_delta_down[i,], nrow = T,  ncol = 1, byrow = TRUE)
# Q_reg_delta_down_day[which(Q_reg_delta_down_day<0)] <- 0
#
# # Q_reg_delta <- matrix(df_soc_q[,alpha], nrow = nrow(df_soc_q)/T,  ncol = T, byrow = TRUE)
# # Q_reg_delta_day <- matrix(Q_reg_delta[i,], nrow = T,  ncol = 1, byrow = TRUE)
# # Q_reg_delta_day_up <- Q_reg_delta_day
# # Q_reg_delta_day_down <- Q_reg_delta_day
# # Q_reg_delta_day_down[which(Q_reg_delta_day_down<0)] <- 0
# # Q_reg_delta_day_up[which(Q_reg_delta_day_up>0)] <- 0
# # Q_reg_delta_day_up


################################ CONSTRAINS #########################################

# constraints <- list(## Model the battery storage
#                     q[1] == Q_start*Q,
#                     q[T+1] == Q_end*Q,
#                     diff(q) == -(bc*eff+bd/eff)*dt+c_reg_share*(Q_reg_delta_down_day-Q_reg_delta_up_day),
#                     -q <= 0,#-Q_min*Q,
#                     q <= Q,
#                     bc <= 0,
#                     -bd <= 0,
#                     b==bc+bd,
#
#                     ## Linearize binary and variable product for bd <= b_ch_dis*(C-c_reg_share*C_reg_up),
#                     z >= 0,
#                     z <= b_ch_dis,
#                     -C_reg_up*(1-b_ch_dis) <= z - c_reg_share,
#                     z - c_reg_share <= 0,
#                     bd <= b_ch_dis*(C)-z*C_reg_up,
#
#                     ## Linearize binary and variable product for bc >= (b_ch_dis-1)*(C-c_reg_share*C_reg_down),
#                     z1 >= 0,
#                     z1 <= b_ch_dis,
#                     -(1-b_ch_dis) <= z1 - c_reg_share,
#                     z1 - c_reg_share <= 0,
#                     bc >= b_ch_dis*(C)-z1*C_reg_down-C+c_reg_share*C_reg_down,
#
#                     ## Split net load in positive and negative for ToU and FiT charge
#                     nl == nl_alpha_day - b + c_reg_share*(Q_reg_delta_down_day-Q_reg_delta_up_day)/dt,
#                     nl == nl_pos + nl_neg,
#                     nl_pos <= nl_pos_neg*nl_fuse, # limits on net-load power
#                     -nl_neg <= -(nl_pos_neg-1)*nl_fuse, # limits on net-load power fuse size
#                     -nl_pos <= 0,
#                     nl_neg <= 0
# )

#################################### OBJECTIVES ########################################

# obj_profit <- (## ToU
#                 t(ptou)%*%(nl_pos-c_reg_share*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt +
#                ## Battery operation
#                 t(pb)%*%(bd)*dt + t(pb)%*%(-bc)*dt + t(pb)%*%((Q_reg_delta_down_day+Q_reg_delta_up_day)*c_reg_share) +
#                ## FiT
#                 t(pfit)%*%(nl_neg-c_reg_share*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt -
#                ## EFR
#                 t(pfcr)%*%(c_reg_share*C)*dt +
#                 t(ptou)%*%(Q_reg_delta_down_day*c_reg_share) + t(pfit)%*%(-Q_reg_delta_up_day*c_reg_share))
#
# const_fcr <- c(constraints, c(q[1:T] <= Q-c_reg_share*C*dt/2,
#                               q[1:T] >= c_reg_share*C*dt/2))
#
# prob_profit <- Problem(Minimize(obj_profit), const_fcr)
# result <- solve(prob_profit, solver = "GUROBI", verbose = FALSE, feastol = 1e-6)
# result$value
#
# obj_env <- t(pco2_alpha_day)%*%(nl)*dt #+ t(pco2_alpha_day)%*%(-bc+solar_alpha_day)*(1-eff^2)*dt
# prob_env <- Problem(Minimize(obj_env), constraints)
# result <- solve(prob_env, solver = "GUROBI", verbose = FALSE, feastol = 1e-6)
# result$value
#
#
# obj_seff <- t(e)%*%(nl_pos-nl_neg)*dt
# prob_seff <- Problem(Minimize(obj_seff), constraints)
# result <- solve(prob_seff, solver = "GUROBI", verbose = FALSE, feastol = 1e-6)
# result$value
#
#
# obj_blackout <- (t(e)%*%(-q[1:T]))
# prob_blackout <- Problem(Minimize(obj_blackout), constraints)
# result <- solve(prob_blackout, solver = "GUROBI", verbose = FALSE, feastol = 1e-6)
# result$value
#
#
# #problems <- c(prob_profit, prob_env, prob_seff, prob_blackout)
# #constrs <- c(const_fcr, constraints, constraints, constraints)
# objectives <- c(obj_profit, obj_env, obj_seff, obj_blackout)
# dm_preferences <- c(0.468+0.01625, 0.216+0.01625, 0.236+0.01625, 0.015+0.01625)
# #dm_preferences <- c(0.25, 0.25, 0.25, 0.25)
# weights <- diag(1, 4, 4)
# result_data <- matrix(, nrow = 4, ncol = 4)
#
#
# for (i in 1:length(objectives)){
#   prob <- Problem(Minimize(weighted_sum(objectives, weights[i,])), constraints)
#   result <- solve(prob, solver = "GUROBI", verbose = FALSE, feastol = 1e-6)
#   print(result$value)
#   ### Fill in the objectives values
#   result_data[i,1] <- ( ## ToU
#                         t(ptou)%*%(result$getValue(nl_pos)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt +
#                         ## Battery operation
#                         t(pb)%*%(result$getValue(bd))*dt + t(pb)%*%(-result$getValue(bc))*dt +
#                         t(pb)%*%((Q_reg_delta_down_day+Q_reg_delta_up_day)*result$getValue(c_reg_share)) +
#                         ## FiT
#                         t(pfit)%*%(result$getValue(nl_neg)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt -
#                         ## EFR
#                         t(pfcr)%*%(result$getValue(c_reg_share)*C)*dt +
#                         t(ptou)%*%(Q_reg_delta_down_day*result$getValue(c_reg_share)) + t(pfit)%*%(-Q_reg_delta_up_day*result$getValue(c_reg_share)))
#
#   result_data[i,2] <- t(pco2_alpha_day)%*%(result$getValue(nl))*dt
#   result_data[i,3] <- t(e)%*%(result$getValue(nl_pos)-result$getValue(nl_neg))*dt
#   result_data[i,4] <- t(e)%*%(-result$getValue(q)[1:T])
# }
#
#
# norm_objectives <- objectives
# for (i in 1:length(objectives)){
#   norm_objectives[[i]] <- get_normalize_obj(result_data, objectives, dm_preferences, i)
# }
#
# prob <- Problem(Minimize(weighted_sum(norm_objectives, c(1.,1.,1.,1.))), const_fcr) #constraints)
# result <- solve(prob, solver = "GUROBI", verbose = FALSE, feastol = 1e-6)
# result$status
# result$solve_time
#
# gresult_data <- matrix(, nrow = 4, ncol = 1)
# for (i in 1:length(gresult_data)){
#   gresult_data[1] <- (## ToU
#                       t(ptou)%*%(result$getValue(nl_pos)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt +
#                       ## Battery operation
#                       t(pb)%*%(result$getValue(bd))*dt + t(pb)%*%(-result$getValue(bc))*dt +
#                       t(pb)%*%((Q_reg_delta_down_day+Q_reg_delta_up_day)*result$getValue(c_reg_share)) +
#                       ## FiT
#                       t(pfit)%*%(result$getValue(nl_neg)-result$getValue(c_reg_share)*((Q_reg_delta_down_day-Q_reg_delta_up_day)/dt))*dt -
#                       ## EFR
#                       t(pfcr)%*%(result$getValue(c_reg_share)*C)*dt +
#                       t(ptou)%*%(Q_reg_delta_down_day*result$getValue(c_reg_share)) + t(pfit)%*%(-Q_reg_delta_up_day*result$getValue(c_reg_share)))
#
#   gresult_data[2] <- t(pco2_alpha_day)%*%(result$getValue(nl))*dt
#   gresult_data[3] <- t(e)%*%(result$getValue(nl_pos)-result$getValue(nl_neg))*dt
#   gresult_data[4] <- t(e)%*%(-result$getValue(q)[1:T])
# }
#
# gresult_data
# result_data
################################### RESULTS #####################################

# result$value
# result$status
# result$getValue(b)
# result$getValue(bc)
# result$getValue(bd)
# result$getValue(c_reg_share)
# result$getValue(z1)
# result$getValue(z)
# result$getValue(b_ch_dis)
# result$getValue(nl)
# result$getValue(nl_pos)
# result$getValue(nl_neg)
# result$getValue(q)
#
# plot(nl_alpha_day, type="l", lwd=2, col="red", asp=2)
# lines(0.5*result$getValue(q),type="b")
# lines(result$getValue(-b), type="l", lwd=2, col="blue", asp=25)
# lines(result$getValue(nl), col="darkred",lwd=2)
# lines(result$getValue(nl_pos), col='darkblue', lwd=3)
# lines(result$getValue(nl_neg), col='lightblue', lwd=3)
# lines(m_real[i,], col='darkgreen', lwd=2)
# lines(m_real[i,]-result$getValue(b), lwd=2)
# lines(4*pco2_alpha_day/max(pco2_alpha_day))

#################################### END #########################################



####### Uncertainty budget

# qd = Variable(T)
# zd = Variable(T)
# y = Variable(T)
#gg <- sample(c(rep(0, 24), rep(1, 24)))
# qd >= 0,
# zd >= 0,
# y >=1,
#constraints <- c(constraints,  p[l] + b[l] == (mu[l]*scaler + qd[l] + gg[l]*zd[l]))
#constraints <- c(constraints,  qd[l] + zd[l] >= mu[l]*scaler*y[l])
# for(l in 2:T){
#   constraints <- c(constraints, q[l]==q[l-1]-(bc[l-1]*eff+bd[l-1]/eff)*dt)# + Q_reg_delta[l]) #c_reg_share
# }
## Garch parameters

# get_robustness_factor <- function(dccrl, W, T, y_fit, r_type) {
#   ## Create Y vector of model fit
#   Y <- list()
#   ## iterate through forecast models
#   for (h in 1:sum(lengths(dccrl@mforecast))){
#     Y_h <- list()
#     z2 <- dim(dccrl@mforecast[[h]]@model[["H"]])[3]
#     z1 <- z2-y_fit+1
#     ## iterate through the fit window
#     for (i in z1:z2){
#       Hi <- dccrl@mforecast[[h]]@model[["H"]][,,i] # matrix of covariances at step i
#       res_i <- dccrl@mforecast[[h]]@model[["residuals"]][i,] # vector of residuals at step i
#       res_i <- matrix(res_i, byrow = TRUE, nrow = T, ncol = 1)
#       Yi <- norm(chol(solve(Hi)) %*% (res_i), type=r_type)
#       Y_h <- append(Y_h, Yi)
#     }
#
#     ## Append the first vector Yi for t+1 to the list
#     Y <- append(Y, list(Y_h[sort.list(unlist(Y_h))]))
#
#     ## Append the vectors for H+t+2:H+t+refit-1
#     ## -1 because we use the new model after refit
#     for (z in 1:(dccrl@model[["out.sample"]][h]-1)){
#       ## remove the first value from the vector
#       Y_h = Y_h[2:length(Y_h)]
#       ## Calculate Yi
#       Hi <- dccrl@mforecast[[h]]@mforecast[["H"]][[z]][,,1] # matrix of covariances at step z
#       mui <- matrix(dccrl@mforecast[[h]]@mforecast[["mu"]][,,z], byrow = TRUE, nrow = T, ncol = 1)
#       id <- z2+z
#
#       idx <- dccrl@model[["rollind"]][[h]][id]
#       print(idx)
#       xi <- matrix(rX[idx,], byrow = TRUE, nrow = T, ncol = 1)
#       Yi = norm(chol(solve(Hi))%*% (xi - mui), type=r_type)
#       Y_h[[length(Y_h)+1]] <- Yi
#
#       Y[[length(Y)+1]] <- Y_h[sort.list(unlist(Y_h))]
#       #Y <- append(Y, list((Y_h[sort.list(unlist(Y_h))])))
#     }
#   }
#   return(t(do.call("cbind", Y)))
# }


# get_mu <- function(dccrl, T) {
#   mx_list <- list()
#   for (h in 1:sum(lengths(dccrl@mforecast))){
#     nrows <- dccrl@model[["out.sample"]][h]
#     mx_list[[length(mx_list)+1]] <- matrix(dccrl@mforecast[[h]]@mforecast$mu, byrow = TRUE, nrow = nrows, ncol = T)
#   }
#   #print(do.call(rbind, mx_list))
#   return(do.call(rbind, mx_list))
# }
#
#
# get_cov <- function(dccrl, T) {
#   "
#   Get aggregated predicted covariance matrix
#   "
#   mx_list <- list()
#   for (h in 1:sum(lengths(dccrl@mforecast))){
#     nrows <- dccrl@model[["out.sample"]][h]
#     mx_list <- append(mx_list, dccrl@mforecast[[h]]@mforecast$H)
#   }
#   #print(do.call(rbind, mx_list))
#   #return(do.call(rbind, mx_list))
#   return(mx_list)
# }

#
# normalize <- function(x) {
#   return((x- min(x)) /(max(x)-min(x)))
# }
#
# denormalize <- function(x,minval,maxval) {
#   x*(maxval-minval) + minval
# }
#
# ddnorm <- as.data.frame(lapply(df[6],normalize))
# ddnorm
# minvec <- sapply(as.data.frame(m),min)
# minvec
# maxvec <- sapply(as.data.frame(m),max)
# maxvec
# dd_n <- as.data.frame(Map(denormalize,as.data.frame(mu_matrix),minvec,maxvec))
# dd_n



# predict_garch <- function(rX, horizon, out_of_sample, T, window, refit){
#   ##----------------------------------------------------------------
#   "Multivariate DCC GARCH forecasting
#    based on http://eclr.humanities.manchester.ac.uk/index.php/R_GARCH"
#   ##----------------------------------------------------------------
#
#   ## specify univariate model
#   uspec.n = multispec(replicate(T, ugarchspec(mean.model = list(armaOrder = c(1,0)))))
#   # specify the DCC model
#   spec1 = dccspec(uspec = uspec.n, dccOrder = c(1, 1), distribution = c("mvnorm", "mvt", "mvlaplace"))
#   dccrl <- dccroll(spec1, data=rX, n.ahead = horizon, forecast.length = out_of_sample,
#                    refit.every = refit, n.start = NULL,
#                    refit.window = c("moving"), window.size = window_size,
#                    solver = "solnp", solver.control = list(),
#                    fit.control = list(eval.se = FALSE))
#   return(dccrl)
# }
#
# predict_quantiles <- function(df, horizon, out_of_sample, T, window_size, refit, N_scenarios){
#   for (colname in colnames(df)){
#     print(colname)
#     "Create dataframe for column"
#     m <- matrix(df[[colname]], byrow = TRUE, nrow = nrow(df)/T, ncol = T)
#     rX = as.data.frame(m)
#     rX['t'] <- seq(as.Date("2013-01-01 00:00:00"), by = "day", length.out = nrow(df)/T)
#     rX = xts(rX[, 1:T], order.by=as.Date(rX$t, format = "%Y-%m-%d", tz = "UTC"))
#
#     "Predict using GARCH"
#     dccrl <- predict_garch(rX, horizon, out_of_sample, T, window_size, refit)
#     mu_matrix <- get_mu(dccrl, T)
#     H_list <- get_cov(dccrl, T)
#     m_real <- m[(window_size+1):dim(m)[1],]
#     print(rmse(c(m_real)-c(mu_matrix)))
#
#     "Create quantiles"
#     quantiles <- seq(from=0.05, to=0.95, by=0.1)
#     m_quantiles <- matrix(, nrow = length(mu_matrix), ncol = length(quantiles)+1)
#     m_quantiles[,length(quantiles)+1] <- matrix(t(mu_matrix), byrow=TRUE, nrow=length(mu_matrix), ncol=1)
#     for (i in 1:dim(m_real)[1]){
#       U <- (chol((H_list[[i]][,,1])))
#       mu <- matrix(c(rep.int(0,T)), nrow = T, ncol = 1, byrow = TRUE)
#       X <- rmvnc(N_scenarios, mu_matrix[i,], U)
#       for (a in 1:length(quantiles)){
#         m_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a],
#                                                         mean=colMeans(X),
#                                                         sd=c(apply(X, 2, sd)))
#
#         # m_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a],
#         #                                                 mean=mu_matrix[i,],
#         #                                                 sd=diag(chol(H_list[[i]][,,1])))
#       }
#     }
#     "Create and save dataframe"
#     df_q <- as.data.frame(m_quantiles)
#     colnames(df_q) <- append(quantiles, 'mean')
#     write.csv(df_q, sprintf('./forecasts/%s_forecast.csv', colname),row.names = FALSE)
#   }
#   return(df_q)
# }



# for (i in 41){
#   X <- diag(x = rep(1, len=T), ncol = T, nrow = T)
#   for (a in seq(from=0.025, to=0.95, by=0.1)){
#     G_a <- Y_matrix[i, round(y_fit*(a))][[1]]
#     #print(a)
#     zn <- c()
#     zp <- c()
#     for(l in 1:T){
#       zn[l] <- mu_matrix[i,l] - G_a*norm(t(solve(chol(solve(H_list[[i]][,,1])))) %*% (X[l,]), type='i')
#       zp[l] <- mu_matrix[i,l] + G_a*norm(t(solve(chol(solve(H_list[[i]][,,1])))) %*% (X[l,]), type='i')
#     }
#     lines((zn), lwd=2)
#     lines((zp), lwd=2)
#   }
#   lines(m_real[i,], col='red',lwd=3)
#   lines(mu_matrix[i,], col='blue',lwd=3)
# }



# for (i in colnames(df)[1]){
#   print(i)
#   "Create dataframe for column i"
#   m <- matrix(df[[i]], byrow = TRUE, nrow = nrow(df)/T, ncol = T)
#   scaler <- max(m)#colMaxs(rX)
#   #m <- m/scaler
#
#   t <- seq(as.Date("2013-01-01 00:00:00"), by = "day", length.out = nrow(df)/T)
#   rX = as.data.frame(m)
#   #rX <- as.data.frame(lapply(rX,normalize))
#   rX['t'] <- t
#   rX = xts(rX[, 1:T], order.by=as.Date(rX$t, format = "%Y-%m-%d", tz = "UTC"))
#
#   ##----------------------------------------------------------------
#   "Multivariate DCC GARCH forecasting
#    based on http://eclr.humanities.manchester.ac.uk/index.php/R_GARCH"
#   ##----------------------------------------------------------------
#
#   ## specify univariate model
#
#   # ugs <- ugarchspec(mean.model = list(include.mean=TRUE),
#   # variance.model = list(garchOrder = c(1,1), model = 'sGARCH'),
#   # distribution.model = 'norm')
#   # uspec.n = multispec(replicate(T, ugs))
#
#   uspec.n = multispec(replicate(T, ugarchspec(mean.model = list(armaOrder = c(1,0)))))
#
#   # specify the DCC model
#   spec1 = dccspec(uspec = uspec.n, dccOrder = c(1, 1), distribution = 'mvnorm')
#   dccrl <- dccroll(spec1, data=rX, n.ahead = 1, forecast.length = out_of_sample,
#                    refit.every = refit, n.start = NULL,
#                    refit.window = c("moving"), window.size = window,
#                    solver = "solnp", solver.control = list(),
#                    fit.control = list(eval.se = FALSE))#, stationarity = FALSE))
#   #Y_matrix <- get_robustness_factor(dccrl, window, T, y_fit, r_type)
#   #res <- get_residuals(dccrl, window, T, y_fit, r_type)
#   mu_matrix <- get_mu(dccrl, T)
#   H_list <- get_cov(dccrl, T)
#   m_real <- m[(window+1):dim(m)[1],]
# }


#q[1] == q[T]-(bc[T]*eff+bd[T]/eff)*dt+Q_reg_delta[T],
#diff(q) == -(bc[1:(T-1)]*eff+bd[1:(T-1)]/eff)*dt+Q_reg_delta[1:(T-1)],# + c_reg_share*Q_reg_delta,

# max_p <- c(0.0, 0.10791766, 0.22484305, 0.3273215 , 0.42466097, 0.52259655, 0.62675316, 0.74546115, 0.8970305 , 1.15691559)
# min_p <- c(0.00614889, 0.10783659, 0.21300555, 0.3202084 , 0.42894743, 0.5392806 , 0.65235873, 0.77227506, 0.91266327, 1.14536361)
# delta_energy <- c(-0.10966215, -0.05396164, -0.02914063, -0.01368981, -0.00331468,
#                   0.00331468,  0.01368981,  0.02914063,  0.05396164,  0.10966215)
# Q_reg_delta <- matrix(c(delta_energy[alpha]*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE)
# C_reg_down = matrix(c(min_p[alpha]*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE)
# C_reg_up = matrix(c(max_p[alpha]*rep.int(1,T)), nrow = T, ncol = 1, byrow = TRUE)
