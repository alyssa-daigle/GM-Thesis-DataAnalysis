# ------------------------------
# Function: ssbyvar
# ------------------------------
# Computes the fraction of total variance explained by a categorical factor
# response      : numeric vector of observations
# category.vec  : factor or character vector indicating groups
# Returns       : numeric fraction (0–1)
# ------------------------------
ssbyvar <- function(response, category.vec) {
  # Mean per category
  means <- tapply(response, category.vec, mean, na.rm = TRUE)

  # Residual sum of squares (difference from category mean)
  ssresid <- sum(sapply(sort(unique(category.vec)), function(z) {
    sum(
      (response[category.vec == z] - means[names(means) == z])^2,
      na.rm = TRUE
    )
  }))

  # Total sum of squares (difference from overall mean)
  sstot <- sum((response - mean(response, na.rm = TRUE))^2, na.rm = TRUE)

  # Fraction of variance explained by the factor
  (sstot - ssresid) / sstot
}
