get_normalize_obj <- function(result_data, objectives, dm_preferences, N){
  max_obj <- max(result_data[,N])
  min_obj <- min(result_data[,N])
  return(dm_preferences[N]*((objectives[[N]]-min_obj)/(max_obj-min_obj)))
}