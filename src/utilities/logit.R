# Logit transformation to use for percentages
logit <- function(value) {
  if (!(all(na.omit(value) >= 0) && all(na.omit(value) <= 1)))
    stop("values must be between 0 and 1.")
  log( value / (1 - value) )
}

# Inverse of the logit transformation
invlogit <- function(value) {
  1 / (1 + exp(-value))
}
