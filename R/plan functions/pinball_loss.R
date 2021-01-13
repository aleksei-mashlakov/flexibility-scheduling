
pinLoss <- function(y, mu, qu, add = TRUE){
  
  # Recursive call for multiple quantiles
  if( length(qu) > 1 ){  
    n <- length( qu )
    l <- sapply(1:n,
                function(ii){
                  return( pinLoss(y, mu[ , ii], qu[ii], add = add) )
                })
    
    if( is.matrix(l) ){ colnames(l) <- qu } else { names(l) <- qu }
    
    return( l )
  }
  
  tau <- 1 - qu
  d <- y - mu
  l <- d * 0
  
  l[d < 0] <- - tau * d[ d < 0 ]
  l[d > 0] <- - (tau-1) * d[ d > 0 ]
  
  if( add ){ l <- sum(l) }
  
  return( l )
}