library(ggplot2)
# Read data
df <- read.csv(here::here("Simulation/sim_data.csv"))

MSE_DiRD <- function(y, x, time_var, c, t0){
  # Prepare regression
  reg_data <- data.frame(y, x, time_var)
  
  # Treatment start indicator
  reg_data$time_indc <- ifelse(reg_data$time_var >= t0, 1, 0)
  
  # Treatment group indicator
  reg_data$group_indc <- ifelse(reg_data$x >= c, 1, 0)
  
  # Interaction between time and treated
  reg_data$dird <- reg_data$time_indc * reg_data$group_indc
  
  # Check if the data is valid for regression
  if (nrow(reg_data) == 0 || sum(!is.na(reg_data$y)) == 0) {
    print("Insufficient data")
  }
  else{
    # Regression model
    model<- lm(y ~ x + group_indc*(1+x) + time_indc*(1+x) + dird*(1+x), data=reg_data)
    mse <- mean(resid(model)^2)
  }
  return(mse)
}

# Function to find the next h value given a current h value
find_next_h <- function(current_h, h_values) {
  # Subset to find all values greater than the current h
  possible_values <- h_values[h_values > current_h]
  
  # Check if there are any higher values
  if(length(possible_values) > 0) {
    return(min(possible_values))
  } else {
    # Find values less than the current h and return the maximum of these
    lower_values <- h_values[h_values < current_h]
    return(max(lower_values))  # Return the maximum of lower values
  }
}


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

DiRD_BW <- function(y, x, c, time_var, t0, ID, max_iter=100){
  ## Initialisation
  
  # Preparation
  # Set up h as the 'bandwidth variable'
  data <- data.frame(ID=ID, y=y, x=x, time_var=time_var)
  data$h <- abs(data$x)
  
  # Initialize parameters : 'for' loop
  h_vector <- c()
  mse_vector <- c()
  log_mse_vector <- c()
  initial_lr <- 0.01
  epsilon <- 1e-8
  
  for(i in 1:max_iter){
    
    # Initialisation : h0 and gradient0 to do a stochastic gradient descent
    gradient_details <- MSE_Gradient0(y=y, x=x, h=data$h, data=data, ID=ID, c=c, time_var=time_var, t0=t0)
    h_0 <- gradient_details[1]
    Gradient_0 <- gradient_details[2]
    h_j <- h_0
    Gradient_j <- Gradient_0
    accumulated_sq_grad <- Gradient_j
    
    # While loop is supposed to find local minima    
    # MSE convergence criterion
    while(abs(Gradient_j) > 0){
      
      # AdaGrad update for Gradient_j
      accumulated_sq_grad <- accumulated_sq_grad + Gradient_j^2
      if (accumulated_sq_grad + epsilon > 0) {  # Check to avoid NaNs from sqrt
        adjusted_lr <- initial_lr / sqrt(accumulated_sq_grad + epsilon)
      } else {
        adjusted_lr <- initial_lr  # Fallback if denominator is zero
      }      
      h_j <- abs(h_j - adjusted_lr * Gradient_j)  # Apply update
      
      # Sort data by unique h values to find the closest h value (h_j_plus)
      unique_h_values <- sort(unique(data$h))
      h_j_plus <- find_next_h(h_j, unique_h_values)
      
      if (h_j_plus < h_j){
        break
      }
      
      # Subset by |x| < h and compute MSE 
      sample_j <- subset(data, subset = abs(x) < h_j)
      sample_j_plus <- subset(data, subset = abs(x) < h_j_plus)
      
      if (nrow(sample_j) == 0 || nrow(sample_j_plus) == 0) {
        # Handle empty data scenario
        next  # Skip to the next iteration
      } else {
        # Calculate MSEs for bothsubsets
        mse_j <- MSE_DiRD(y = sample_j$y, x=sample_j$x, 
                          time_var = sample_j$time_var, c=c, t0 = t0)
        mse_j_plus <- MSE_DiRD(y = sample_j_plus$y, x=sample_j_plus$x, 
                               time_var = sample_j_plus$time_var, c=c, t0 = t0)
      }
      
      # Compute Gradient
      Gradient_j <- (mse_j - mse_j_plus) / max((h_j - h_j_plus), 1e-05)
    }
    
    # h_j is the local minimizer for MSE
    h_vector <- c(h_vector, h_j)
    mse_vector <- c(mse_vector, mse_j)
    log_mse_vector <- c(log_mse_vector, log(mse_j + 1))  # Log-transform MSE to stabilize variance
  }

  local_minima <- data.frame(h=h_vector, mse=mse_vector, log_mse=log_mse_vector)
  # Find the h corresponding to the minimal mse
  minimal_mse_h <- local_minima$h[which.min(local_minima$log_mse)]
  
  mse_plot <- ggplot(local_minima, aes(x = h, y = log_mse)) +
    geom_line() +
    geom_point() +
    labs(title = "Evolution of MSE over h values", x = "Bandwidth (h)", y = "MSE") +
    theme_minimal()
  
  return(list( optimal_h = minimal_mse_h, 
              min(local_minima$log_mse),
              plot = mse_plot))
}

xd<- DiRD_BW(y=df$Price_m2, x=df$Distance_to_Border, c=0, time_var=df$Year, t0=3,
             ID=df$PropertyID)
xd
