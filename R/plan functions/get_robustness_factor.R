get_robustness_factor <- function(dccrl, W, k, y_fit, r_type) {
  ## Create Y vector of model fit
  Y <- list()
  ## iterate through forecast models
  for (h in 1:sum(lengths(dccrl@mforecast))){
    Y_h <- list()
    z2 <- dim(dccrl@mforecast[[h]]@model[["H"]])[3]
    z1 <- z2-y_fit+1
    ## iterate through the fit window
    for (i in z1:z2){
      Hi <- dccrl@mforecast[[h]]@model[["H"]][,,i] # matrix of covariances at step i
      res_i <- dccrl@mforecast[[h]]@model[["residuals"]][i,] # vector of residuals at step i
      res_i <- matrix(res_i, byrow = TRUE, nrow = k, ncol = 1)
      Yi <- norm(chol(solve(Hi)) %*% (res_i), type=r_type)
      Y_h <- append(Y_h, Yi)
    }

    ## Append the first vector Yi for t+1 to the list
    Y <- append(Y, list(Y_h[sort.list(unlist(Y_h))]))

    ## Append the vectors for H+t+2:H+t+refit-1
    ## -1 because we use the new model after refit
    for (z in 1:(dccrl@model[["out.sample"]][h]-1)){
      ## remove the first value from the vector
      Y_h = Y_h[2:length(Y_h)]
      ## Calculate Yi
      Hi <- dccrl@mforecast[[h]]@mforecast[["H"]][[z]][,,1] # matrix of covariances at step z
      mui <- matrix(dccrl@mforecast[[h]]@mforecast[["mu"]][,,z], byrow = TRUE, nrow = k, ncol = 1)
      id <- z2+z

      idx <- dccrl@model[["rollind"]][[h]][id]
      print(idx)
      xi <- matrix(rX[idx,], byrow = TRUE, nrow = k, ncol = 1)
      Yi = norm(chol(solve(Hi))%*% (xi - mui), type=r_type)
      Y_h[[length(Y_h)+1]] <- Yi

      Y[[length(Y)+1]] <- Y_h[sort.list(unlist(Y_h))]
      #Y <- append(Y, list((Y_h[sort.list(unlist(Y_h))])))
    }
  }
  return(t(do.call("cbind", Y)))
}