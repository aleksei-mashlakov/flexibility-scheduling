optimizee <- function(I_step) {
  best_Ns <- list()
  best_Nsmses <- list()
  for (l in 1:(I_step-1)){
    G_a <- Y_matrix[l, round(y_fit*(0.975))][[1]]
    N_rmse_best <- Inf
    bestN <- 1
    for (N in 1:y_fit){
      if (N==1){
        mu_N <- (mu_matrix[l,]+G_a*(solve(chol(solve(H_list[[l]][,,1]))))%*%t(t(do.call(cbind, res[l,])[,1])))*scaler
        N_rmse <- rmse(m_real[l,]*scaler-mu_N)
      }
      else{
        mu_N <- (mu_matrix[l,]+G_a*(solve(chol(solve(H_list[[l]][,,1]))))%*%t(t(rowMeans(do.call(cbind, res[l,])[,1:N]))))*scaler
      }
      N_rmse <- rmse(m_real[l,]*scaler-mu_N)
      if (N_rmse < N_rmse_best){
        N_rmse_best <- N_rmse
        bestN <- N
      }
    }
    #print(bestN)
    #Y_h <- append(Y_h, Yi)
    best_Nsmses[[length(best_Nsmses)+1]] <- N_rmse_best
    best_Ns[[length(best_Ns)+1]] <- bestN
    # zn <- c()
    # zp <- c()
    # for(l in 1:T){
    #   zn[l] <- mu[l] - G_a*norm(t(solve(chol(solve(H_list[[i]][,,1])))) %*% (X[l,]), type='I')
    #   zp[l] <- mu[l] + G_a*norm(t(solve(chol(solve(H_list[[i]][,,1])))) %*% (X[l,]), type='I')
    # }
    # lines((zn)*scaler, lwd=2)
    # lines((zp)*scaler, lwd=2)
  }
  print(mean(do.call(cbind, best_Nsmses)))
  print(round(mean(do.call(cbind, best_Ns))))
  return(round(mean(do.call(cbind, best_Ns))))
}