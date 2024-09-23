data <- df

mini_batch_gradient <- function(y, x, h, h_0, data, ID, c, time_var, t0, mini_batch_size=35){

  # Sort data by closeness to h_0
  data_unique <- data[!duplicated(data$PropertyID),]
  unique_h_values <- data_unique$h
  
  mini_batch_size <- min(nrow(data_unique), mini_batch_size)  # Ensures not to exceed data size
  
  mini_batch_h_vector <- mini_batch_select(data=data_unique, h_0=h_0, mini_batch_size=mini_batch_size)
  
  gradient_vector <- c()
  
  for (hi in mini_batch_h_vector) {
    # Get h_plus vector
    h_plus <- find_next_h(current_h = hi, h_values = mini_batch_h_vector)
    
    # Subset by |x| < h and compute MSE (for h_0)
    eval_sample <- subset(data, subset = h < hi)
    # Repeat the process for h_plus
    sample_plus <- subset(data, subset = h < h_plus)
    
    if (nrow(eval_sample) == 0 || nrow(sample_plus) == 0) {
      # Skip to the next iteration if data is empty
      next  
    } else {
      # Calculate MSEs for bothsubsets
      mse_0 <- MSE_DiRD(y = eval_sample$Price_m2, x=eval_sample$Distance_to_Border, 
                        time_var = eval_sample$Year, c=0, t0 = 3)
      mse_plus <- MSE_DiRD(y = sample_plus$Price_m2, x=sample_plus$Distance_to_Border, 
                           time_var = sample_plus$Year, c=0, t0 = 3)
    }
    # Gradient
    MSE_gradient <-  (mse_0-mse_plus)/max((h_0-h_plus), 1e-08)
    gradient_vector <- c(gradient_vector, MSE_gradient)
  }
  return(mean(gradient_vector))
}

mini_batch_gradient(y=data$Price_m2, x=data$Distance_to_Border, h=data$h, h_0=h0
                    data=data, ID= data$PropertyID, c=0, time_var=data$Year, t0=3)
