MSE_Gradient0 <- function(y, x, h, data, ID, c, time_var, t0){
  # get a df with all unique IDs
  data_unique <- data[!duplicated(data$ID),]
  unique_h_values <- data_unique$h
  
  # Randomly pick a random h and the one just above (h_plus)
  random_index <- sample(1:nrow(data_unique), 1)
  h_0 <- data_unique$h[random_index] # rnd h
  h_plus <- find_next_h(current_h = h_0, h_values = unique_h_values) # h just above
  
  # Subset by |x| < h and compute MSE (for h_0)
  eval_sample <- subset(data, subset = h < h_0)
  # Repeat the process for h_plus
  sample_plus <- subset(data, subset = h < h_plus)
  
  
  if (nrow(eval_sample) == 0 || nrow(sample_plus) == 0) {
    # Handle empty data scenario
    # If data is insufficient, re-run the function
    return(MSE_Gradient0(y, x, h, data, ID, c, time_var, t0))  # Skip to the next iteration
  } else {
    # Calculate MSEs for bothsubsets
    mse_0 <- MSE_DiRD(y = eval_sample$y, x=eval_sample$x, 
                      time_var = eval_sample$time_var, c=c, t0 = t0)
    mse_plus <- MSE_DiRD(y = sample_plus$y, x=sample_plus$x, 
                         time_var = sample_plus$time_var, c=c, t0 = t0)
  }
  # Gradient
  MSE_gradient <-  (mse_0-mse_plus)/max((h_0-h_plus), 1e-08)
  return(c(h_0, MSE_gradient))
}