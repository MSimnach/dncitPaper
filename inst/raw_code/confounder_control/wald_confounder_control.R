# Load necessary library
library(MASS)
library(RCIT)

# Simulation parameters
n <- 1000  # Number of samples
p <- 20  # Number of dimensions for X
num_simulations <- 40  # Number of simulations

# Store results
wald_p_values <- numeric(num_simulations)
rcot_p_values <- numeric(num_simulations)

# Function to simulate data and perform Wald test
simulate_study <- function() {
  # Generate confounder Z
  Z <- rnorm(n, mean=0, sd=1)

  # Generate a p-dimensional vector of X as a function of Z plus noise
  eps_X <- mvrnorm(n, mu=rep(0, p), Sigma=diag(p))
  X <- matrix(rep(Z, p), ncol=p) + matrix(rep(Z^2, p), ncol=p) + eps_X

  # Generate the outcome Y with a nonlinear effect of Z and no effect of X
  eps <- rnorm(n, mean=0, sd=1)
  Y <- Z+Z^3 + eps  # Nonlinear function of Z, no effect of X

  # Fit the linear model Y ~ X + Z
  colnames(X) <- paste0("X", 1:p)
  data <- as.data.frame(cbind(Y, X, Z))
  model <- lm(Y ~ . -1, data=data)  # Exclude intercept to match dimension

  wtest <- lmtest::waldtest(model, colnames(X))
  p_value_wald <- wtest$`Pr(>F)`[2]

  # Perform RCOT
  p_value_rcot <- RCIT::RCoT(Y, X, Z)$p
  return(c(p_value_wald, p_value_rcot))
}

# Perform the simulations
for (i in 1:num_simulations) {
  if(i %% 20 == 0){print(i)}
  p_values <- simulate_study()
  wald_p_values[i] <- p_values[1]
  rcot_p_values[i] <- p_values[2]
}

# Analyze the results
significance_level <- 0.05
proportion_significant_wald <- mean(wald_p_values < significance_level)
proportion_significant_rcot <- mean(rcot_p_values < significance_level)

# Output the results
cat("Proportion of significant results:", proportion_significant, "\n")
hist(wald_p_values, breaks=30, main="Histogram of Wald Test P-values", xlab="P-value")
hist(rcot_p_values, breaks=30, main="Histogram of RCOT P-values", xlab="P-value")
abline(v=0.05, col="red", lwd=2)
