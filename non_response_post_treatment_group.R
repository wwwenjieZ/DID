

############################################################

#####***prost-treatment non response

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

data_potential <- data%>%
  mutate(non_response_prob = pnorm(Y, mean = 0, sd = 1))%>%
  mutate(non_response = ifelse(non_response_prob > 0.9, 1, 0)) %>%
  filter(!(time == 1 & C>0 & non_response == 1)) %>%
  select(id, time, D, C, Y)
nrow(data_potential)# 8045

model_results_OLS_AT <- vector("list", n_simulations)

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
      
      # Define treatment and potential outcomes
      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
      D = ifelse(Dtrue , 1, 0) 
      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
                           ifelse(time == 0 & D == 1, U_0_1,
                                  ifelse(time == 1 & D == 0, U_1_0,
                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
      
      # Simulate outcomes Y based on potential outcomes
      Y <- alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected
    })
  
  # Simulate non-response
  data_AT <- data %>%
    mutate(non_response_prob = pnorm(Y, mean = 0, sd = 1))%>%
    mutate(non_response = ifelse(non_response_prob > 0.9, 1, 0)) %>%
    filter(!(time == 1 & C>0 & non_response == 1)) %>%
    select(id, time, D, C, Y,d)
  
  
  data_AT <- data_AT %>% mutate(did = d * time)
  didOLS_AT <- lm(Y ~ d + time + did, data = data_AT)
  model_results_OLS_AT[[i]] <- coef(didOLS_AT)
}

#----results---#
results_matrix_OLS_AT <- do.call(rbind, model_results_OLS_AT)
mean_coeffs_OLS_AT <- apply(results_matrix_OLS_AT, 2, mean)
sd_coeffs_OLS_AT <- apply(results_matrix_OLS_AT , 2, sd)
bias_OLS_AT <- mean_coeffs_OLS_AT - true_coeffs_OLS
rmse_OLS_AT <- sqrt(apply((results_matrix_OLS_AT - true_coeffs_OLS)^2, 2, mean))

#----Print the results---#
print(mean_coeffs_OLS_AT)
print(sd_coeffs_OLS_AT)
print(bias_OLS_AT)
print(rmse_OLS_AT)


####---------------------FE--------------------######

model_results_FE_AT <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2) )
  
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
  
  # Simulate non-response
  data_AT <- data %>%
    mutate(non_response_prob = pnorm(Y, mean = 0, sd = 1))%>%
    mutate(non_response = ifelse(non_response_prob > 0.9, 1, 0)) %>%
    filter(!(time == 1 & C > 0 & non_response == 1)) %>%
    select(id, time, D, C, Y,d)
  
  
  data_AT <- data_AT %>% mutate(did = d * time)
  didFE_AT <-plm(Y ~ d+ time+did, data = data_AT, model = "within")
  model_results_FE_AT[[i]] <- coef(didFE_AT)
}

#----results---#
results_matrix_FE_AT <- do.call(rbind, model_results_FE_AT)
mean_coeffs_FE_AT <- apply(results_matrix_FE_AT, 2, mean)
sd_coeffs_FE_AT <- apply(results_matrix_FE_AT , 2, sd)
bias_FE_AT <- mean_coeffs_FE_AT - true_coeffs_FE
rmse_FE_AT <- sqrt(apply((results_matrix_FE_AT - true_coeffs_FE)^2, 2, mean))

#----Print the results---#
print(mean_coeffs_FE_AT)
print(sd_coeffs_FE_AT)
print(bias_FE_AT)
print(rmse_FE_AT)
