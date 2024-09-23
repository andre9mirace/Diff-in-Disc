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