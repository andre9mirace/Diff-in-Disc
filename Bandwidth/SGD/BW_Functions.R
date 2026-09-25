# Bandwidth selection functions for Diff-in-Disc.
# Documented in BW_Functions.Rmd, which displays each section with knitr::read_chunk().

## ---- Imbens_BW ----
Imbens_BW <- function(y, x, c, kernel="triangular"){
  
  df <- data.frame(y=y, x=x)
  df <- na.omit(df)
  N <- nrow(df)
  
  # Uniform kernel
  if (kernel == "uniform") {
   C_k <- 5.40
  } 
    # Triangular kernel as default
  else {
    C_k <- 3.4375
  }

  # Prepare the df for regression
  df$group_indc <- ifelse(df$x>c, 1, 0) # group indicator for regression
  df$x1 <- (df$x-c)
  df$x2 <- (df$x-c)^2
  df$x3 <- (df$x-c)^3
  
  # Subset into left/right of the cutoff
  sample_plus <- subset(df, x > c)
  sample_minus <- subset(df, x < c)
  N_plus <- nrow(sample_plus)
  N_minus <- nrow(sample_minus)
  
  # First step
  # We first estimate the running variable's density and the limit of the conditional variances at the cutoff 
  
  h1 <- 1.84*sd(x)*N^(-1/5) # Standard Silverman rule of thumb
  silver_sample_plus <- subset(sample_plus, x < c+h1)
  silver_sample_minus <- subset(sample_minus, x > c-h1)
  
  running_density_c <- (nrow(silver_sample_minus)+nrow(silver_sample_plus))/(2*N*h1) # Density of X at the cutoff
  
  # Conditional variance limits (c+ and c-)
  cond_var_c_plus <- var(silver_sample_plus$y)
  cond_var_c_minus <- var(silver_sample_minus$y)
  
  # Second step
  # Estimation of two bandwidths to estimate the second order derivative.

  # Quadratic regression
  reg_m3 <- lm(y ~ group_indc + x1 + x2 + x3, data = df) # maybe adding time_indc for DiRD?
  
  m3_c <- 6*coef(reg_m3)[[5]] # estimate of the 3rd derivative to the regression function 

  
  h2_plus <- 3.56*((cond_var_c_plus)/(running_density_c*(m3_c)^2*N_plus))^(1/7)
  h2_minus <- 3.56*((cond_var_c_minus)/(running_density_c*(m3_c)^2*N_minus))^(1/7)
  
  sample_2_plus <- subset(sample_plus, sample_plus$x < c+h2_plus)
  sample_2_minus <- subset(sample_minus, sample_minus$x > c-h2_minus)
  
  
  # Quadratic regression to estimate 
  quad_reg_plus <- lm(y ~ x1 + x2, data = sample_2_plus)
  quad_reg_minus <- lm(y ~ x1 + x2, data = sample_2_minus)
  
  m2_plus <- 2*coef(quad_reg_plus)[[3]]
  m2_minus <- 2*coef(quad_reg_minus)[[3]]

  # Third step
  r_plus <- 2160*cond_var_c_plus/(nrow(sample_2_plus)*h2_plus^4)
  r_minus <- 2160*cond_var_c_minus/(nrow(sample_2_minus)*h2_minus^4)
  
  h_optimal <- C_k*N^(-1/5)*((cond_var_c_plus+cond_var_c_minus)
                             /(running_density_c*(m2_plus-m2_minus)^2+r_plus+r_minus))^(1/5)
  return(h_optimal)
}

## ---- rdd_k_weights ----
rdd_k_weights <- function(x, c, h, kernel="triangular"){
  # Normalize the distances from the cutoff by dividing by the bandwidth (h).
  u = (x - c) / h
  
	# Uniform kernel
  if (kernel == "uniform") {
   w = (0.5 * (abs(u) <= 1)) / h
  } 
    # Triangular kernel as default
  else {
     w = ((1 - abs(u)) * (abs(u) <= 1)) / h
  }
  return(w)  
}

## ---- MSE_DiRD ----
MSE_DiRD <- function(y, x, time_var, c, t0, h, kernel = "triangular"){
  # Keep observations inside the bandwidth (they all have a positive kernel weight)
  keep <- !is.na(y) & !is.na(x) & abs(x - c) < h
  y <- y[keep]
  x <- x[keep]
  time_var <- time_var[keep]

  # Calculate weights
  w <- rdd_k_weights(x=x, c=c, h=h, kernel = kernel)

  # Treatment start indicator, treatment group indicator and their interaction
  time_indc <- as.numeric(time_var >= t0)
  group_indc <- as.numeric(x >= c)
  dird <- time_indc * group_indc

  # Check if every group x period cell has enough data for the regression
  if (any(tabulate(1 + group_indc + 2*time_indc, nbins = 4) < 3)) {
    return(NA)
  }

  # Same design as lm(y ~ x + group_indc*(1+x) + time_indc*(1+x) + dird*(1+x), weights = w),
  # fitted with .lm.fit because this function is called thousands of times by DiRD_BW
  X <- cbind(1, x, group_indc, time_indc, dird, x*group_indc, x*time_indc, x*dird)
  fit <- .lm.fit(X * sqrt(w), y * sqrt(w))

  # Unweighted residuals, as resid() returns for a weighted lm
  mse <- mean((fit$residuals / sqrt(w))^2)
  return(mse)
}

## ---- MSE_penalty ----
MSE_penalty <- function(N, n_h, penalty_type="default"){
  
  if (penalty_type == "default") {
    penalty <- (N / n_h)
    
    # Log
  } else if (penalty_type == "log") {
    penalty <- log(N / n_h + 1)
    
    # Sqrt penalty
  } else if (penalty_type == "sqrt") {
    penalty <- sqrt(N / n_h)
    
    # Cube root pen
  } else if (penalty_type == "cbrt") {  # Cube root
    penalty <- 0.01*(N / (n_h))^(1/3)
    
    # Exp penalty
  } else if (penalty_type == "exp" | penalty_type == "exponential") {
    alpha <- 0.01  
    penalty <- exp(alpha * N / n_h) - 1
  } else if (penalty_type %in% c("sigmoid", "sig")) {
    alpha <- 0.05  
    beta <- 10     
    penalty <- 1 / (1 + exp(-alpha * (N / n_h - beta))) # sigmoid function
  } 
  
  return(penalty)
}

## ---- penalized_MSE ----
penalized_MSE <- function(data, h, c, t0, N, kernel = "triangular", penalty_type = "default"){
  sample_h <- data[abs(data$x - c) < h, ]
  if (nrow(sample_h) == 0) {
    return(NA)
  }
  penalty <- MSE_penalty(N = N, n_h = nrow(sample_h), penalty_type = penalty_type)
  mse <- MSE_DiRD(y = sample_h$y, x = sample_h$x, time_var = sample_h$time_var,
                  c = c, t0 = t0, h = h, kernel = kernel)
  return(penalty * mse)
}

## ---- MSE_gradient ----
# Central finite difference of the penalized MSE with a step proportional to h
MSE_gradient <- function(data, h, c, t0, N, kernel = "triangular", penalty_type = "default",
                         step = 0.02){
  delta <- step * h
  mse_up <- penalized_MSE(data = data, h = h + delta, c = c, t0 = t0, N = N,
                          kernel = kernel, penalty_type = penalty_type)
  mse_down <- penalized_MSE(data = data, h = h - delta, c = c, t0 = t0, N = N,
                            kernel = kernel, penalty_type = penalty_type)
  return((mse_up - mse_down) / (2 * delta))
}

## ---- MSE_Gradient0 ----
MSE_Gradient0 <- function(data, h_candidates, c, t0, N,
                          kernel="triangular", penalty_type="default"){

  for (attempt in 1:100) {
    # Randomly pick a starting bandwidth
    h_0 <- h_candidates[sample.int(length(h_candidates), 1)]

    # Penalized MSE and its gradient at h_0
    mse_0 <- penalized_MSE(data = data, h = h_0, c = c, t0 = t0, N = N,
                           kernel = kernel, penalty_type = penalty_type)
    Gradient_0 <- MSE_gradient(data = data, h = h_0, c = c, t0 = t0, N = N,
                               kernel = kernel, penalty_type = penalty_type)

    # If data is insufficient, draw another h_0
    if (!is.na(mse_0) & !is.na(Gradient_0)) {
      return(c(h_0, mse_0, Gradient_0))
    }
  }
  return(c(NA, NA, NA))
}

## ---- DiRD_BW ----
DiRD_BW <- function(y, x, c, time_var, t0, ID, 
                    kernel = "triangular", penalty_type = "default",
                    max_iter=10, max_epochs=50){
  ## Initialization

  # Preparation
  data <- data.frame(ID=ID, y=y, x=x, time_var=time_var)
  data$h <- abs(data$x - c)
  N <- nrow(data)
  h_values <- sort(unique(data$h)) # every admissible bandwidth
  # Initialize parameters : 'for' loop
  h0_vector <- c()
  h_vector <- c()
  mse_vector <- c()
  initial_lr <- 1
  epsilon <- 1e-8

  # h intervals based on quantiles of h
  quantile_division <- c()
  if(N >= 4000){
    quantile_division <- c(0, 0.05, 0.15, 0.25, 0.40, 0.55, 0.75, 1)
  } else if (N < 300) {
    quantile_division <- c(0, 0.45, 1)
  } else {
      quantile_division <- seq(from = 0, to = 1, by = 0.1)
    }
  h_quantiles <- quantile(data$h, probs = quantile_division)

  # Loop through each quantile interval
  for (q in 1:(length(h_quantiles) - 1)) {
    h_min <- h_quantiles[q]
    h_max <- h_quantiles[q + 1]

    # The interval only restricts where h_0 is drawn: the regressions always use
    # the full data with |x - c| < h (overwriting `data` here emptied every
    # interval after the first one)
    h_candidates <- h_values[h_values >= h_min & h_values < h_max]
    if (length(h_candidates) == 0) next  # Skip if no data in the current interval

    # AdaGrad base learning rate: a tenth of the interval width, in the unit of x
    eta <- initial_lr * (h_max - h_min) / 10

    for(i in 1:max_epochs){

      # Initialization : h0 and gradient0 to do a stochastic gradient descent
      gradient_details <- MSE_Gradient0(data = data, h_candidates = h_candidates,
                                        c = c, t0 = t0, N = N,
                                        kernel = kernel, penalty_type = penalty_type)
      if (anyNA(gradient_details)) next
      h_0 <- gradient_details[1]

      # Prepare for iterations (mse_j always matches h_j)
      h_j <- h_0
      mse_j <- gradient_details[2]
      Gradient_j <- gradient_details[3]
      accumulated_sq_grad <- 0

      # Best bandwidth visited during this descent
      h_best <- h_j
      mse_best <- mse_j

      unchanged_iter_count <- 0  # Counter for unchanged h_j

      for(iter in (1:max_iter)){

        # AdaGrad update: eta over the root of the sum of all squared gradients
        accumulated_sq_grad <- accumulated_sq_grad + Gradient_j^2
        adjusted_lr <- eta / sqrt(accumulated_sq_grad + epsilon)

        # SGD Update
        h_j_new <- h_j - adjusted_lr * Gradient_j

        # Stay within the range of admissible bandwidths
        if (h_j_new <= min(h_values) | h_j_new >= max(h_values)) break

        # Penalized MSE and its gradient at the new bandwidth
        mse_j_new <- penalized_MSE(data = data, h = h_j_new, c = c, t0 = t0, N = N,
                                   kernel = kernel, penalty_type = penalty_type)
        Gradient_j <- MSE_gradient(data = data, h = h_j_new, c = c, t0 = t0, N = N,
                                   kernel = kernel, penalty_type = penalty_type)

        # Handle insufficient data scenario
        if (is.na(mse_j_new) | is.na(Gradient_j)) {
          break
        }

        if (mse_j_new < mse_best) {
          h_best <- h_j_new
          mse_best <- mse_j_new
        }

        # Threshold
        if (abs(mse_j_new-mse_j) < 1/log(N)){
          break
        }

        # Check if h_j remains unchanged
        if (h_j == h_j_new) {
          unchanged_iter_count <- unchanged_iter_count + 1
        } else {
          unchanged_iter_count <- 0  # Reset counter if h_j changes
        }

        if (unchanged_iter_count >= 2) {
          break  # Break if h_j doesn't change for 2 iterations
        }

        h_j <- h_j_new  # Update h_j to the new value
        mse_j <- mse_j_new
      }
      # store starting point, local minimizer and minimums for MSE
      h0_vector <- c(h0_vector, h_0)
      h_vector <- c(h_vector, h_best)
      mse_vector <- c(mse_vector, mse_best)

    }
  }
  local_minima <- data.frame(h_0=h0_vector, h=h_vector, mse=mse_vector)
  local_minima <- na.omit(local_minima)

  minimal_mse_h <- local_minima$h[which.min(local_minima$mse)]

  return(list(optimal_h = minimal_mse_h, optimal_mse = min(local_minima$mse),
              local_minima = local_minima))
}
