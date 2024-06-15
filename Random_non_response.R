

############################################################

#####***Random non-response

###########################################################

##################***OLS########
rm(list = ls())
library(dplyr)
library(plm)
library(dplyr)

set.seed(123)
n <- 5000 
alpha <- 1
gamma <- -1
beta1 <- 2
beta2 <- -2
theta <- 1
delta <- 5 
pR2 <- 0.7
n_simulations=100
model_results_OLS_random <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2)  
  )
  
  data <- data |>  
    within({  
      C <- rnorm(n * 2, mean = 0, sd = 1)    
      U_0_0 <- rnorm(n * 2, mean = 0, sd = 1)  
      U_0_1 <- rnorm(n * 2, mean = 0, sd = 1)  
      U_1_0 <- rnorm(n * 2, mean = 0, sd = 1)  
      U_1_1 <- rnorm(n * 2, mean = 0, sd = 1)  
      E <- rnorm(n * 2, mean = 0, sd = 1)  
      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
      D = ifelse(Dtrue , 1, 0) 
      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
                           ifelse(time == 0 & D == 1, U_0_1,
                                  ifelse(time == 1 & D == 0, U_1_0,
                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
      
      # Simulate outcomes Y based on potential outcomes
      Y <- alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected
    })
  
  data_pR2_sim <- data  
  D0 <- which(data_pR2_sim$D == 0)  
  D1 <- which(data_pR2_sim$D == 1)  
  
  #Randomly remove some data to form a new non-response data set
  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
  data_pR2_sim <- data_pR2_sim[-rows_to_remove, ] 
  
  #  OLS regression
  data_pR2_sim <- data_pR2_sim %>% mutate(did = d * time)
  didOLS_random <- lm(Y ~ d + time + did, data = data_pR2_sim)
  model_results_OLS_random[[i]] <- coef(didOLS_random)
}

#----results---#
results_matrix_OLS_random <- do.call(rbind, model_results_OLS_random)
mean_coeffs_OLS_random <- apply(results_matrix_OLS_random, 2, mean)
sd_coeffs_OLS_random <- apply(results_matrix_OLS_random , 2, sd)
bias_OLS_random <- mean_coeffs_OLS_random - true_coeffs_OLS
rmse_OLS_random <- sqrt(apply((results_matrix_OLS_random - true_coeffs_FE)^2, 2, mean))

#----Print the results---#
print(mean_coeffs_OLS_random)
print(sd_coeffs_OLS_random)
print(bias_OLS_random)
print(rmse_OLS_random)

####---------------------FE--------------------######

model_results_FE_random <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2)  
  )
  
  data <- data |>  
    within({  
      C <- rnorm(n * 2, mean = 0, sd = 1)    
      U_0_0 <- rnorm(n * 2, mean = 0, sd = 1)  
      U_0_1 <- rnorm(n * 2, mean = 0, sd = 1)  
      U_1_0 <- rnorm(n * 2, mean = 0, sd = 1)  
      U_1_1 <- rnorm(n * 2, mean = 0, sd = 1)  
      E <- rnorm(n * 2, mean = 0, sd = 1)  
      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
      D = ifelse(Dtrue , 1, 0) 
      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
                           ifelse(time == 0 & D == 1, U_0_1,
                                  ifelse(time == 1 & D == 0, U_1_0,
                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
      
      # Simulate outcomes Y based on potential outcomes
      Y <- alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected
    })
  
  data_pR2_sim <- data  
  D0 <- which(data_pR2_sim$D == 0)  
  D1 <- which(data_pR2_sim$D == 1)  
  
  #Randomly remove some data to form a new non-response data set
  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
  data_pR2_sim <- data_pR2_sim[-rows_to_remove, ] 
  
  #  FE regression
  data_pR2_sim <- data_pR2_sim %>% mutate(did = d * time)
  didFE_random <- plm(Y ~ d+ time+did, data = data_pR2_sim, model = "within")
  model_results_FE_random[[i]] <- coef(didFE_random)
}

#----results---#
results_matrix_FE_random <- do.call(rbind, model_results_FE_random)
mean_coeffs_FE_random <- apply(results_matrix_FE_random, 2, mean)
sd_coeffs_FE_random <- apply(results_matrix_FE_random , 2, sd)
bias_FE_random <- mean_coeffs_FE_random - true_coeffs_FE
rmse_FE_random <- sqrt(apply((results_matrix_FE_random - true_coeffs_FE)^2, 2, mean))

#----Print the results---#
print(mean_coeffs_FE_random)
print(sd_coeffs_FE_random)
print(bias_FE_random)
print(rmse_FE_random)
