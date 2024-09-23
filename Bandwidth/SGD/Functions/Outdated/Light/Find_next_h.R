# Function to find the next h value given a current h value
find_next_h <- function(current_h, h_values) {
  # Ensure h_values is not empty
  if (length(h_values) == 0) {
    stop("h_values cannot be empty")
  }
  # Subset to find all values greater than the current h
  possible_values <- h_values[h_values > current_h]
  
  # If there are any higher values, return the minimum
  if (length(possible_values) > 0) {
    return(min(possible_values))
  } else {
    # Find values less than the current h and return the maximum of these
    lower_values <- h_values[h_values < current_h]
    if (length(lower_values) > 0) {
      return(max(lower_values))
      
      # No higher or lower values available = stop
    } else {
      stop("No higher or lower values found for the given h")
    }
  }
}