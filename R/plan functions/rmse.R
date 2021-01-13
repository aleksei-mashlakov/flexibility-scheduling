rmse <- function(x) {
    sqrt(mean(x^2, na.rm = TRUE))
}