# Data generating process of the Paris / Issy-les-Moulineaux simulation
# (Simulation/data_generation.Rmd) written as a function, so it can be drawn many times.
#
# Price_m2 = 10000 + 500 Balcony - 1.7 Distance_to_EiffelTower + 200 Metro_access
#            + 50 Local_Trend + tau Treatment + noise
#
# design = "sharp"  : local-trend zone with a sharp edge at |x| = border_range (the thesis design).
#                     Inside the zone the trend starts at 100 and grows 4% a year; outside it starts
#                     at 200 and grows 8% (Paris side) or 6% (Issy side).
# design = "smooth" : the same trends, but the move from the border zone to the outer zones is
#                     gradual, s(x) = 1 - exp(-(x / border_range)^2), so the trend gap between
#                     both sides grows smoothly (quadratically near the border) with the distance.
#
# thesis_rng = TRUE reproduces Simulation/sim_data.csv with seed = 123. data_generation.Rmd calls
# set.seed() before X, Y and the balcony draw, which makes the three perfectly dependent
# (X = 0.8 (Y + 2000), Balcony = 1 iff Y > 1500). thesis_rng = FALSE draws them independently,
# as the report describes.
simulate_sim_data <- function(seed, design = c("sharp", "smooth"), n_properties = 10000,
                              n_years = 6, t0 = 3, tau = 350, border_range = 500,
                              noise_level = 50, thesis_rng = FALSE){
  design <- match.arg(design)
  N <- n_properties * n_years

  # Coordinates and balcony (one draw per property)
  set.seed(seed)
  X_axis <- runif(n_properties, min = 0, max = 4000)
  if (thesis_rng) set.seed(seed)
  Y_axis <- runif(n_properties, min = -2000, max = 3000)
  if (thesis_rng) set.seed(seed)
  balcony <- sample(c(0, 1), size = n_properties, replace = TRUE, prob = c(0.7, 0.3))

  data <- data.frame(PropertyID = rep(1:n_properties, each = n_years),
                     Year = rep(1:n_years, n_properties),
                     X_axis = rep(X_axis, each = n_years),
                     Y_axis = rep(Y_axis, each = n_years))
  data$Region <- ifelse(data$Y_axis > 0, "Paris", "Issy")
  data$Distance_to_Border <- data$Y_axis
  data$Treatment <- ifelse(data$Region == "Paris" & data$Year >= t0, 1, 0)
  data$Balcony_indc <- rep(balcony, each = n_years)
  data$Distance_to_EiffelTower <- sqrt(data$X_axis^2 + (data$Y_axis - 3400)^2)
  data$Metro_access <- ifelse(data$Region == "Paris", 1, 0)

  # Local trend: level in year 1 and yearly growth rate
  x <- data$Distance_to_Border
  outer_growth <- ifelse(x > 0, 0.08, 0.06)
  if (design == "sharp") {
    s <- as.numeric(abs(x) > border_range)
  } else {
    s <- 1 - exp(-(x / border_range)^2)
  }
  initial_level <- 100 + 100 * s
  growth <- 0.04 + (outer_growth - 0.04) * s
  data$Local_Trend <- initial_level * (1 + growth)^(data$Year - 1)

  # Noise, twice as strong after year 3 (drawn in the same order as data_generation.Rmd)
  data$asymmetric_noise <- ifelse(data$Year > 3, rnorm(N, 0, noise_level * 2),
                                  rnorm(N, 0, noise_level))

  data$Price_m2 <- 10000 + 500 * data$Balcony_indc - 1.7 * data$Distance_to_EiffelTower +
    200 * data$Metro_access + 50 * data$Local_Trend + tau * data$Treatment +
    data$asymmetric_noise
  data
}
