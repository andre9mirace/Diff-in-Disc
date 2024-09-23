DiRD_BW <- function(y, x, c, time_var, t0, ID, max_iter=10, max_epochs=200){
  ## Initialisation
  
  # Preparation
  data <- data.frame(ID=ID, y=y, x=x, time_var=time_var)
  data$h <- abs(data$x)
  N <- nrow(data)
  # Initialize parameters : 'for' loop
  h_vector <- c()
  mse_vector <- c()
  initial_lr <- 0.01
  epsilon <- 1e-8
  
  # h intervals based on quantiles of h
  quantile_division <- c()
  if(N >= 4000){
    quantile_division <- c(0, 0.05, 0.15, 0.25, 0.4, 0.55, 0.75, 1)
  } else if (N < 300) {
    quantile_division <- c(0, 0.3, 1)
  } else {
      quantile_division <- seq(from = 0, to = 1, by = 0.1)
    }
  h_quantiles <- quantile(data$h, probs = quantile_division)

  # Adapter par quantiles. Voir, une division plus adaptée a petits jeux de données et/ou grand db. Eviter les cas de petits subsets
  # Plusieurs h_0 par intervalle -> plus de chances de converger dans les petites valeurs. 
  # boucle de h_0 en general
  # Loop through each quantile interval
  for (q in 1:(length(h_quantiles) - 1)) {
    h_min <- h_quantiles[q]
    h_max <- h_quantiles[q + 1]
    
    # Subset data within the current h quantile interval
    data <- subset(data, h >= h_min & h < h_max)
    if (nrow(data) == 0) next  # Skip if no data in the current interval
    
    # Get h values in the current interval
    current_h_values <- sort(unique(data$h))
    
    # Initialise empty mse_j
    mse_j <- NA
    mse_j_plus <- NA
    
    for(i in 1:max_epochs){
      
      # Initialisation : h0 and gradient0 to do a stochastic gradient descent
      gradient_details <- MSE_Gradient0(y=y, x=x, h=data$h, data=data, ID=ID, c=c, time_var=time_var, t0=t0)
      h_0 <- gradient_details[1]
      Gradient_j <- gradient_details[2]
      
      # Prepare for iterations
      h_j <- h_0
      iter_count <- 0
      accumulated_sq_grad <- Gradient_j^2
      
      # Dynamic threshold
      all_gradients <- c()
      # After calculating Gradient_j, add it to all_gradients
      all_gradients <- c(all_gradients, abs(Gradient_j))
      edf <- ecdf(all_gradients)
      
      # Calculate dynamic threshold based on the standard deviation of recorded gradients
      if (length(all_gradients) > 1) {
        dynamic_threshold <- quantile(edf, 0.8 * max_epochs - i / max_epochs) # As per TA algorithm (Gilli & al.)
      } else {
        dynamic_threshold <- 10  # Default threshold for the first epoch
      }
      
      while(abs(Gradient_j) > dynamic_threshold){
        
        # AdaGrad update for Gradient_j
        accumulated_sq_grad <- accumulated_sq_grad + Gradient_j^2
        adjusted_lr <- initial_lr / sqrt(accumulated_sq_grad + epsilon)
        
        h_j <- h_j - adjusted_lr * Gradient_j

        # Sort data by unique h values to find the closest h value (h_j_plus)
        unique_h_values <- data[!duplicated(data$ID),]
        unique_h_values <- sort(unique(data$h))
        
        h_j_plus <- find_next_h(h_j, unique_h_values)
        
        # Handle hj max scenario
        if (h_j_plus < h_j){
          next
        }
        
        # Subset by |x| < h and compute MSE 
        sample_j <- subset(data, subset = abs(x-c) < h_j)
        sample_j_plus <- subset(data, subset = abs(x-c) < h_j_plus)
       
        # Handle empty data scenario
        if (nrow(sample_j) == 0 | nrow(sample_j_plus) == 0) {
          mse_j <- NA
          
        } else {
          # Calculate MSEs for bothsubsets
          mse_j <- MSE_DiRD(y = sample_j$y, x=sample_j$x, 
                            time_var = sample_j$time_var, c=c, t0 = t0)
          mse_j_plus <- MSE_DiRD(y = sample_j_plus$y, x=sample_j_plus$x, 
                                 time_var = sample_j_plus$time_var, c=c, t0 = t0)
  
          # Compute Gradient
          Gradient_j <- (mse_j - mse_j_plus) / max((h_j - h_j_plus), 1e-05)
        }
        
        # Break when you get to the max_iterations
        if (iter_count >= max_iter){
          break
        }
  
      }
      # store local minimizer and minimums for MSE
      h_vector <- c(h_vector, h_j)
      mse_vector <- c(mse_vector, mse_j)
    }
  }
  local_minima <- data.frame(h=h_vector, mse=mse_vector)
  local_minima <- na.omit(local_minima)
  
  minimal_mse_h <- local_minima$h[which.min(local_minima$mse)]
  
  return(list
         (optimal_h = minimal_mse_h, min(local_minima$mse))
         )
}

DiRD_BW(y=df$Price_m2, x=df$Distance_to_Border, c=0, time_var=df$Year, t0=3, ID=df$PropertyID)



bw_vect <- c()

for (val in (1:100)){
  bw <- DiRD_BW(y=df$Price_m2, x=df$Distance_to_Border, c=0, time_var=df$Year, t0=3, ID=df$PropertyID)
  bw_vect <- c(bw_vect, bw)
}
