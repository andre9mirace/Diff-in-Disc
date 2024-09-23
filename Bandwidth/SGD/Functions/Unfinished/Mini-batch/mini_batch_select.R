mini_batch_select <- function(data, h_0, mini_batch_size) {
  data$distance <- abs(data$h - h_0)  # Compute absolute differences
  data_sorted <- data[order(data$distance), ]  # Sort by distance
  mini_batch <- head(data_sorted, mini_batch_size) # Select the closest 'b' obs
  return(mini_batch$h) # return the vector of bandwidths
}
