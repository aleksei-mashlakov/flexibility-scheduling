predict_garch <- function(rX, horizon, out_of_sample, T, window_size, refit){
  ##----------------------------------------------------------------
  "Multivariate DCC GARCH forecasting
   based on http://eclr.humanities.manchester.ac.uk/index.php/R_GARCH"
  ##----------------------------------------------------------------
  
  ## specify univariate model
  uspec.n = rugarch::multispec(replicate(T, rugarch::ugarchspec(mean.model = list(armaOrder = c(1,0)))))
  # specify the DCC model
  spec1 = rmgarch::dccspec(uspec = uspec.n, dccOrder = c(1, 1), distribution = c("mvnorm", "mvt", "mvlaplace"))
  dccrl <- rmgarch::dccroll(spec1, data=rX, n.ahead = horizon, forecast.length = out_of_sample, 
                             refit.every = refit, n.start = NULL,
                             refit.window = c("moving"), window.size = window_size,
                             solver = "solnp", solver.control = list(),
                             fit.control = list(eval.se = FALSE))
  return(dccrl)
}