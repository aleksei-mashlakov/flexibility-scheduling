get_mu <- function(dccrl, k) {
  mx_list <- list()
  for (h in 1:sum(lengths(dccrl@mforecast))){
    nrows <- dccrl@model[["out.sample"]][h]
    mx_list[[length(mx_list)+1]] <- matrix(dccrl@mforecast[[h]]@mforecast$mu, byrow = TRUE, nrow = nrows, ncol = k)
  }
  #print(do.call(rbind, mx_list))
  return(do.call(rbind, mx_list))
}

