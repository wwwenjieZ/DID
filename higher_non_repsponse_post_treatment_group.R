

############################################################

#####***prost-treatment non response unbalance

###########################################################

##################***OLS########
rm(list = ls())
library(dplyr)
library(plm)
library(dplyr)


###########################################
#Data Simulation
######################################


n <- 5000 
alpha <- 1
gamma <- -1
beta1 <- 2
beta2 <- -2
theta <- 1
delta <- 5 

set.seed(123)
data <- data.frame(
  id = rep(1:(n), each = 2),   
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
    D_true  <- ifelse(time == 1, 1, 0)
    D<- D_true 
    U_selected <- ifelse(time == 0 & D == 0, U_0_0,
                         ifelse(time == 0 & D == 1, U_0_1,
                                ifelse(time == 1 & D == 0, U_1_0,
                                       ifelse(time == 1 & D == 1, U_1_1, NA))))
    
    # Simulate outcomes Y based on potential outcomes
    Y <- alpha + beta1 * time + beta2 * D_true + theta * (D_true * d) + C * D_true + U_selected
  })

# Check the data structure 
str(data)
head(data)
sum(data$D == 0)
sum(data$D == 1)

##################################################

#"true” estimator

##################################################
data$time <- as.numeric(data$time)
#OLS estimate
data <- data %>%mutate(post = ifelse(time == 1, 1, 0),did = d * post)
didreg_OLS <- lm(Y ~d+ time+ did, data = data)

#FE  estimate
data <- pdata.frame(data, index=c("id", "time"))
didreg_FE <-plm(Y ~ d+ time+did, data = data, model = "within")
true_coeffs_OLS<-coef(didreg_OLS)
true_coeffs_FE<-coef(didreg_FE)


###############################################

#non-response simulation checking

###############################################
data <- data %>%
  mutate(Y_standardized = (Y - mean(Y)) / sd(Y))
data_potential <- data%>%
  mutate(non_response_prob = pnorm(Y_standardized, mean = 0, sd = 1))%>%
  mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
  filter(!(time == 1 & non_response == 1)) %>%
  select(id, time, D, C, Y,non_response_prob)

nrow(data_potential)
ids_in_both_periods <- intersect(
  data_potential %>% filter(time == 0) %>% select(id) %>% unlist(),
  data_potential %>% filter(time == 1) %>% select(id) %>% unlist()
)
data_potential<-data_potential %>% filter(id %in% ids_in_both_periods)
nrow(data_potential)


###############################################

#unbalance data set

###############################################
n_simulations=1000
set.seed(123)
model_results_OLS_unbal <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2)
  )
  
  data <- data %>%
    mutate(
      C = rnorm(n * 2, mean = 0, sd = 1),
      U_0_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_0_1 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_1 = rnorm(n * 2, mean = 0, sd = 1),
      E = rnorm(n * 2, mean = 0, sd = 1)
    ) %>%
    mutate(
      Dtrue = time > 0 & (gamma + beta1 + beta2 + theta + delta * C + E) > 0,
      D = ifelse(Dtrue, 1, 0),
      U_selected = case_when(
        time == 0 & D == 0 ~ U_0_0,
        time == 0 & D == 1 ~ U_0_1,
        time == 1 & D == 0 ~ U_1_0,
        time == 1 & D == 1 ~ U_1_1
      ),
      Y = alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected,
      did = d * time
    )
  data <- data %>%
    mutate(Y_standardized = (Y - mean(Y)) / sd(Y))
  data <- data %>%
    mutate(non_response_prob = pnorm(Y_standardized , mean = 0, sd = 1)) %>%
    mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
    filter(!(time == 1  & non_response == 1)) %>%
    select(id, time, D, C, Y, d, did)
  
  didOLS_unbal <- lm(Y ~ d + time + did, data = data)
  model_results_OLS_unbal[[i]] <- coef(didOLS_unbal)
}

# Convert results to a data frame
results_df_OLS_unbal <- do.call(rbind, model_results_OLS_unbal) %>%
  as.data.frame()

#calculate bias, std, and RMSE
calculate_metrics_unbal <- function(estimates, true_coeffs_OLS) {
  bias <- mean(estimates - true_coeffs_OLS)
  std <- sd(estimates)
  rmse <- sqrt(mean((estimates - true_coeffs_OLS)^2))
  return(c(bias = bias, std = std, rmse = rmse))
}

metrics_OLS_unbal <- lapply(names(results_df_OLS_unbal), function(param) {
  calculate_metrics_unbal(results_df_OLS_unbal[[param]], true_coeffs_OLS[param])
})

# Convert to data frame and print result
metrics_df_OLS_unbal <- do.call(rbind, metrics_OLS_unbal)
rownames(metrics_df_OLS_unbal) <- names(results_df_OLS_unbal)
metrics_df_OLS_unbal

###############################################

#balance data set

###############################################
set.seed(123)
model_results_FE <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2)
  )
  
  data <- data %>%
    mutate(
      C = rnorm(n * 2, mean = 0, sd = 1),
      U_0_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_0_1 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_1 = rnorm(n * 2, mean = 0, sd = 1),
      E = rnorm(n * 2, mean = 0, sd = 1)
    ) %>%
    mutate(
      Dtrue = time > 0 & (gamma + beta1 + beta2 + theta + delta * C + E) > 0,
      D = ifelse(Dtrue, 1, 0),
      U_selected = case_when(
        time == 0 & D == 0 ~ U_0_0,
        time == 0 & D == 1 ~ U_0_1,
        time == 1 & D == 0 ~ U_1_0,
        time == 1 & D == 1 ~ U_1_1
      ),
      Y = alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected,
      did = d * time
    )
  data <- data %>%
    mutate(Y_standardized = (Y - mean(Y)) / sd(Y))
  data <- data %>%
    mutate(non_response_prob = pnorm(Y_standardized, mean = 0, sd = 1)) %>%
    mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
    filter(!(time == 1  & non_response == 1)) %>%
    select(id, time, D, C, Y, d, did)
  # Create a balanced sample
  ids_in_both_periods <- intersect(
    data %>% filter(time == 0) %>% select(id) %>% unlist(),
    data %>% filter(time == 1) %>% select(id) %>% unlist()
  )
  
  data<-data %>% filter(id %in% ids_in_both_periods)
  data <- pdata.frame(data, index = c("id", "time"))
  didFE <- plm(Y ~ d + time + did, data = data, model = "within")
  model_results_FE[[i]] <- coef(didFE)
}

results_FE <- do.call(rbind, model_results_FE) %>%
  as.data.frame()

# calculate bias, std, and RMSE
calculate_metrics_FE <- function(estimates, true_coeffs_FE) {
  bias <- mean(estimates - true_coeffs_FE)
  std <- sd(estimates)
  rmse <- sqrt(mean((estimates - true_coeffs_FE)^2))
  return(c(bias = bias, std = std, rmse = rmse))
}


metrics_FE <- lapply(names(results_FE), function(param) {
  calculate_metrics_FE(results_FE[[param]], true_coeffs_FE[param])
})


set.seed(123)
model_results_OLS <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2)
  )
  
  data <- data %>%
    mutate(
      C = rnorm(n * 2, mean = 0, sd = 1),
      U_0_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_0_1 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_1 = rnorm(n * 2, mean = 0, sd = 1),
      E = rnorm(n * 2, mean = 0, sd = 1)
    ) %>%
    mutate(
      Dtrue = time > 0 & (gamma + beta1 + beta2 + theta + delta * C + E) > 0,
      D = ifelse(Dtrue, 1, 0),
      U_selected = case_when(
        time == 0 & D == 0 ~ U_0_0,
        time == 0 & D == 1 ~ U_0_1,
        time == 1 & D == 0 ~ U_1_0,
        time == 1 & D == 1 ~ U_1_1
      ),
      Y = alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected,
      did = d * time
    )
  data <- data %>%
    mutate(Y_standardized = (Y - mean(Y)) / sd(Y))
  data <- data %>%
    mutate(non_response_prob = pnorm(Y_standardized , mean = 0, sd = 1)) %>%
    mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
    filter(!(time == 1  & non_response == 1)) %>%
    select(id, time, D, C, Y, d, did)
  # Create a balanced sample
  ids_in_both_periods <- intersect(
    data %>% filter(time == 0) %>% select(id) %>% unlist(),
    data %>% filter(time == 1) %>% select(id) %>% unlist()
  )
  
  data<-data %>% filter(id %in% ids_in_both_periods)
  didOLS <- lm(Y ~ d + time + did, data = data)
  model_results_OLS[[i]] <- coef(didOLS)
}

# Convert results to a data frame
results_OLS <- do.call(rbind, model_results_OLS) %>%
  as.data.frame()

calculate_metrics_OLS <- function(estimates, true_coeffs_OLS) {
  bias <- mean(estimates - true_coeffs_OLS)
  std <- sd(estimates)
  rmse <- sqrt(mean((estimates - true_coeffs_OLS)^2))
  return(c(bias = bias, std = std, rmse = rmse))
}

metrics_OLS <- lapply(names(results_OLS), function(param) {
  calculate_metrics_OLS(results_OLS[[param]], true_coeffs_OLS[param])
})


# print result
metrics_FE <- do.call(rbind, metrics_FE)
rownames(metrics_FE) <- names(results_FE)
metrics_FE

metrics_df_OLS <- do.call(rbind, metrics_OLS)
rownames(metrics_df_OLS) <- names(results_OLS)
metrics_df_OLS




#model_results_OLS_AT <- vector("list", n_simulations)
#
#for (i in 1:n_simulations) {
#  data <- data.frame(
#    id = rep(1:n, each = 2),   
#    time = rep(0:1, times = n),  
#    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2)  
#  )
#  
#  data <- data |>  
#    within({  
#      C <- rnorm(n * 2, mean = 0, sd = 1)    
#      U_0_0 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_0_1 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_1_0 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_1_1 <- rnorm(n * 2, mean = 0, sd = 1)  
#      E <- rnorm(n * 2, mean = 0, sd = 1)  
#      
#      # Define treatment and potential outcomes
#      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
#      D = ifelse(Dtrue , 1, 0) 
#      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
#                           ifelse(time == 0 & D == 1, U_0_1,
#                                  ifelse(time == 1 & D == 0, U_1_0,
#                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
#      
#      # Simulate outcomes Y based on potential outcomes
#      Y <- alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected
#    })
#  
#  # Simulate non-response
#  data_AT <- data %>%
#    mutate(non_response_prob = pnorm(Y, mean = 0, sd = 1))%>%
#    mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
#    filter(!(time == 1 & non_response == 1)) %>%
#    select(id, time, D, C, Y,d)
#  
#  
#  data_AT <- data_AT %>% mutate(did = d * time)
#  didOLS_AT <- lm(Y ~ d + time + did, data = data_AT)
#  model_results_OLS_AT[[i]] <- coef(didOLS_AT)
#}
#
##----results---#
#results_matrix_OLS_AT <- do.call(rbind, model_results_OLS_AT)
#mean_coeffs_OLS_AT <- apply(results_matrix_OLS_AT, 2, mean)
#sd_coeffs_OLS_AT <- apply(results_matrix_OLS_AT , 2, sd)
#bias_OLS_AT <- mean_coeffs_OLS_AT - true_coeffs_OLS
#rmse_OLS_AT <- sqrt(apply((results_matrix_OLS_AT - true_coeffs_OLS)^2, 2, mean))
#
##----Print the results---#
#print(mean_coeffs_OLS_AT)
#print(sd_coeffs_OLS_AT)
#print(bias_OLS_AT)
#print(rmse_OLS_AT)
#
#
#######################
###balance
#####v####################
###---------------------FE--------------------######
#n_simulations=100
#model_results_FE_AT_bal <- vector("list", n_simulations)
#for (i in 1:n_simulations) {
#  data <- data.frame(
#    id = rep(1:n, each = 2),   
#    time = rep(0:1, times = n),  
#    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2) )
#  
#  data <- data |>  
#    within({  
#      C <- rnorm(n * 2, mean = 0, sd = 1)    
#      U_0_0 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_0_1 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_1_0 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_1_1 <- rnorm(n * 2, mean = 0, sd = 1)  
#      E <- rnorm(n * 2, mean = 0, sd = 1)  
#      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
#      D = ifelse(Dtrue , 1, 0) 
#      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
#                           ifelse(time == 0 & D == 1, U_0_1,
#                                  ifelse(time == 1 & D == 0, U_1_0,
#                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
#      
#      # Simulate outcomes Y based on potential outcomes
#      Y <- alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected
#    })
#  
#  # Simulate non-response
#  data_AT <- data %>%
#    mutate(non_response_prob = pnorm(Y, mean = 0, sd = 1))%>%
#    mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
#    filter(!(time == 1  & non_response == 1)) %>%
#    select(id, time, D, C, Y,d)
#  # Create a balanced sample
#  ids_in_both_periods <- intersect(
#    data_AT %>% filter(time == 0) %>% select(id) %>% unlist(),
#    data_AT %>% filter(time == 1) %>% select(id) %>% unlist()
#  )
#  balanced_dataAT <-data_AT %>% filter(id %in% ids_in_both_periods)
#  
#  # FE regression on balanced data
#  balanced_dataAT <- balanced_dataAT %>% mutate(did = d * time)
#  balanced_dataAT_P <- pdata.frame(balanced_dataAT, index=c("id", "time"))
#  didFE_AT_bal <- plm(Y ~ d + time + did, data = balanced_dataAT_P, model="within")
#  model_results_FE_AT_bal[[i]] <- coef(didFE_AT_bal)
#}
#
##----results---#
#results_matrix_FE_AT_bal <- do.call(rbind, model_results_FE_AT_bal)
#mean_coeffs_FE_AT_bal <- apply(results_matrix_FE_AT_bal, 2, mean)
#sd_coeffs_FE_AT_bal <- apply(results_matrix_FE_AT_bal , 2, sd)
#bias_FE_AT_bal <- mean_coeffs_FE_AT_bal - true_coeffs_FE
#rmse_FE_AT_bal <- sqrt(apply((results_matrix_FE_AT_bal - true_coeffs_FE)^2, 2, mean))
#
##----Print the results---#
#print(mean_coeffs_FE_AT_bal)
#print(sd_coeffs_FE_AT_bal)
#print(bias_FE_AT_bal)
#print(rmse_FE_AT_bal)
#
#
#model_results_OLS_AT_bal <- vector("list", n_simulations)
#for (i in 1:n_simulations) {
#  data <- data.frame(
#    id = rep(1:n, each = 2),   
#    time = rep(0:1, times = n),  
#    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2) )
#  
#  data <- data |>  
#    within({  
#      C <- rnorm(n * 2, mean = 0, sd = 1)    
#      U_0_0 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_0_1 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_1_0 <- rnorm(n * 2, mean = 0, sd = 1)  
#      U_1_1 <- rnorm(n * 2, mean = 0, sd = 1)  
#      E <- rnorm(n * 2, mean = 0, sd = 1)  
#      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
#      D = ifelse(Dtrue , 1, 0) 
#      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
#                           ifelse(time == 0 & D == 1, U_0_1,
#                                  ifelse(time == 1 & D == 0, U_1_0,
#                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
#      
#      # Simulate outcomes Y based on potential outcomes
#      Y <- alpha + beta1 * time + beta2 * D + theta * (D * d) + C * D + U_selected
#    })
#  
#  # Simulate non-response
#  data_AT <- data %>%
#    mutate(non_response_prob = pnorm(Y, mean = 0, sd = 1))%>%
#    mutate(non_response = ifelse(non_response_prob > 0.1, 1, 0)) %>%
#    filter(!(time == 1  & non_response == 1)) %>%
#    select(id, time, D, C, Y,d)
#  # Create a balanced sample
#  ids_in_both_periods <- intersect(
#    data_AT %>% filter(time == 0) %>% select(id) %>% unlist(),
#    data_AT %>% filter(time == 1) %>% select(id) %>% unlist()
#  )
#  balanced_dataAT <-data_AT %>% filter(id %in% ids_in_both_periods)
#  
#  # FE regression on balanced data
#  balanced_dataAT <- balanced_dataAT %>% mutate(did = d * time)
#  didOLS_AT_bal <- lm(Y ~ d + time + did, data = balanced_dataAT)
#  model_results_OLS_AT_bal[[i]] <- coef(didOLS_AT_bal)
#}
#
##----results---#
#results_matrix_OLS_AT_bal <- do.call(rbind, model_results_OLS_AT_bal)
#mean_coeffs_OLS_AT_bal <- apply(results_matrix_OLS_AT_bal, 2, mean)
#sd_coeffs_OLS_AT_bal <- apply(results_matrix_OLS_AT_bal , 2, sd)
#bias_OLS_AT_bal <- mean_coeffs_OLS_AT_bal - true_coeffs_OLS
#rmse_OLS_AT_bal <- sqrt(apply((results_matrix_OLS_AT_bal - true_coeffs_OLS)^2, 2, mean))
#
##----Print the results---#
#print(mean_coeffs_OLS_AT_bal)
#print(sd_coeffs_OLS_AT_bal)
#print(bias_OLS_AT_bal)
#print(rmse_OLS_AT_bal)
#
#
#