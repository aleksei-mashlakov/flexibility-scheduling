get_cov <- function(dccrl, k) {
  "
  Get aggregated predicted covariance matrix
  "
  mx_list <- list()
  for (h in 1:sum(lengths(dccrl@mforecast))){
    nrows <- dccrl@model[["out.sample"]][h]
    mx_list <- append(mx_list, dccrl@mforecast[[h]]@mforecast$H)
  }
  #print(do.call(rbind, mx_list))
  #return(do.call(rbind, mx_list))
  return(mx_list)
}