library(MASS)
library(partykit)
library(model4you)
library(palmtree)


#set.seed(123)

n <- 300      # Number of patients
m <- 30       # Total number of variables
p <- 2        # Predictive factors
q <- 2        # Prognostic factors
Delta_beta <- 1
gamma <- c(1,1)

# Generate treatment indicator
X_A <- rbinom(n, size = 1, prob = 0.5)
X_A

# Create correlation matrix Σ
Sigma <- matrix(0.2, m, m)
diag(Sigma) <- 1
Sigma

# Simulate multivariate normal Z
Z <- mvrnorm(n, mu = rep(0, m), Sigma = Sigma)

# Assign variable names
colnames(Z) <- paste0("Z", 1:m)
Z

# Split into categories
Z_predictive <- Z[, 1:p]              # Z1, Z2
Z_prognostic <- Z[, (p+1):(p+q)]      # Z3, Z4
Z_noise <- Z[, (p+q+1):m]             # Z5 to Z30


# Define individual predictive variables
Z1 <- Z_predictive[, 1]
Z2 <- Z_predictive[, 2]

# Define X_F as prognostic variables
X_F <- Z_prognostic


# Define β(Z) based on tree logic
beta_Z <- ifelse(Z1 <= 0, -0.375,
                 ifelse(Z1 > 0 & Z2 <= 0, -0.375 + Delta_beta,
                        -0.375 + (2*Delta_beta)))  # else Z1 > 0 & Z2 > 0


# Simulate error term
U <- rnorm(n, mean = 0, sd = sqrt(1.5))

# Compute outcome Y
Y <- X_A * beta_Z + X_F %*% gamma + U


###############################
# Combine everything into a data frame
#sim_data <- data.frame(
#  Y = Y,
#  a = factor(X_A),
#  Z1 = Z1,
#  Z2 = Z2,
#  Z3 = Z_prognostic[, 1],
#  Z4 = Z_prognostic[, 2]
#)
###############################
# Fit PALM tree - 1
#palm_model <- palmtree(
#  Y ~ a | Z3 + Z4 | Z1 + Z2,
#  data = sim_data
#)
###############################
# Fit palm tree - 2
# Create the formula components
#global_vars <- c("Z3", "Z4")  # Prognostic
#partition_vars <- paste0("Z", 1:30)  # All possible for partitioning

#partition_formula <- as.formula(
#  paste("Y ~ a |", paste(global_vars, collapse = " + "),
#        "|", paste(partition_vars, collapse = " + "))
#)

#palm_model <- palmtree(
#  formula = partition_formula,
#  data = sim_data
#)
###############################

# Reattach all Z1 to Z30, but keep Y and treatment
sim_data <- data.frame(
  Y = as.numeric(Y),  # ensure it's numeric
  a = factor(X_A),      # treatment
  Z                     # Z1 to Z30
)

global_vars <- c("Z3", "Z4")              # Prognostic
partition_vars <- paste0("Z", 1:30)       # Predictive + Noise

library(partykit)
library(model4you)
library(palmtree)

# Define the formula
partition_formula <- as.formula(
  paste("Y ~ a |", paste(global_vars, collapse = " + "), "|", paste(partition_vars, collapse = " + "))
)

# Fit the model
palm_model <- palmtree(
  formula = partition_formula,
  data = sim_data
)

# Inspect results
print(palm_model)
plot(palm_model)



# Coefficients
coef(palm_model, model = "tree")
coef(palm_model, model = "palm")
coef(palm_model, model = "all")

#########################################
library(mclust)  # for adjustedRandIndex

# True group labels based on beta_Z
true_group <- as.factor(beta_Z)

# Extract terminal node (leaf) IDs as predicted clusters from the tree
predicted_group <- as.factor(predict(palm_model, type = "node"))

# Compute Adjusted Rand Index
adjustedRandIndex(true_group, predicted_group)





