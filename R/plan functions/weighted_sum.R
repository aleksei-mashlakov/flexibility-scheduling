weighted_sum <- function(objective, weights) {
  num_objs <- length(objective)
  weighted_objs <- lapply(1:num_objs, function(i) { objective[[i]] * weights[i] })
  Reduce("+", weighted_objs)
}