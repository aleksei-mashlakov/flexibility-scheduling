

predict_quantiles <- function(df, horizon, out_of_sample, T, window_size, refit, N_scenarios){
  " Predict quantiles using GARCH and random sampling
    Inputs: df - class dataframe
  "
  numCores <- parallel::detectCores()
  print(numCores)
  doParallel::registerDoParallel(numCores)

  #numCores <- parallel::detectCores()
  #print(numCores)
  #cl <- parallel::makeCluster(numCores)
  #doParallel::registerDoParallel(cl)

  foreach(col_idx=1:dim(df)[2], .export=c("predict_garch", "get_mu", "get_cov", "rmse"),
                                .packages=c("foreach","doParallel","xts","LaplacesDemon","rmgarch","rugarch")) %dopar% {
    # export=ls(envir=globalenv())
    colname <- colnames(df)[col_idx]
    print(colname)

    "Create dataframe for column"
    m <- matrix(df[[colname]], byrow = TRUE, nrow = nrow(df)/T, ncol = T)
    rX = as.data.frame(m)
    rX['t'] <- seq(as.Date("2013-01-01 00:00:00"), by = "day", length.out = nrow(df)/T)
    rX = xts::xts(rX[, 1:T], order.by=as.Date(rX$t, format = "%Y-%m-%d", tz = "UTC"))

    "Predict using GARCH"
    dccrl <- predict_garch(rX, horizon, out_of_sample, T, window_size, refit)
    mu_matrix <- get_mu(dccrl, T)
    H_list <- get_cov(dccrl, T)
    m_real <- m[(window_size+1):dim(m)[1],]
    print(rmse(c(m_real)-c(mu_matrix)))

    "Create quantiles"
    #quantiles <- seq(from=0.05, to=0.95, by=0.1)
    quantiles <- seq(from=0.95, to=0.05, by=-0.05)
    m_quantiles <- matrix(, nrow = length(mu_matrix), ncol = length(quantiles)+1)
    m_quantiles[,length(quantiles)+1] <- matrix(t(mu_matrix), byrow=TRUE, nrow=length(mu_matrix), ncol=1)
    for (i in 1:dim(m_real)[1]){
      U <- (chol((H_list[[i]][,,1])))
      #mu <- matrix(c(rep.int(0,T)), nrow = T, ncol = 1, byrow = TRUE)
      X <- LaplacesDemon::rmvnc(N_scenarios, mu_matrix[i,], U)
      for (a in 1:length(quantiles)){
        m_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a],
                                                        mean=colMeans(X),
                                                        sd=c(apply(X, 2, sd)))

        # m_quantiles[((i-1)*T+1):((i-1)*T+T),a] <- qnorm(p=quantiles[a],
        #                                                 mean=mu_matrix[i,],
        #                                                 sd=diag(chol(H_list[[i]][,,1])))
      }
    }
    "Create and save dataframe"
    df_q <- as.data.frame(m_quantiles)
    colnames(df_q) <- append(quantiles, 'mean')
    write.csv(df_q, sprintf('./forecasts/%s_forecast.csv', colname),row.names = FALSE)
  }
  #return(df_q)
  # When you're done, clean up the cluster
  doParallel::stopImplicitCluster()
}
