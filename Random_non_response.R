

############################################################

#####***Random non-response-Un-balance

###########################################################

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
    Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected,
    did = d * time
  )

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

#unbalance data set simulation 

##############################################

#pR2=0.5
pR2=0.2
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
      Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected,
      did = d * time
    )
  
   D0 <- which(data$D == 0)  
   D1 <- which(data$D == 1)  
   
   #Randomly remove some data to form a new non-response data set
   n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
   n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
   rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
   rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
   rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
   data <- data[-rows_to_remove, ] 
  
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
  calculate_metrics_unbal (results_df_OLS_unbal[[param]], true_coeffs_OLS[param])
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
      Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected,
      did = d * time
    )
  
  data <- data  
  D0 <- which(data$D == 0)  
  D1 <- which(data$D == 1)  
  
  #Randomly remove some data to form a new non-response data set
  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
  data <- data[-rows_to_remove, ] 
  
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
calculate_metrics <- function(estimates, true_coeffs_FE) {
  bias <- mean(estimates - true_coeffs_FE)
  std <- sd(estimates)
  rmse <- sqrt(mean((estimates - true_coeffs_FE)^2))
  return(c(bias = bias, std = std, rmse = rmse))
}


metrics <- lapply(names(results_FE), function(param) {
  calculate_metrics(results_FE[[param]], true_coeffs_FE[param])
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
      Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected,
      did = d * time
    )
  
  data <- data  
  D0 <- which(data$D == 0)  
  D1 <- which(data$D == 1)  
  
  #Randomly remove some data to form a new non-response data set
  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
  data <- data[-rows_to_remove, ] 
  
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
metrics_FE <- do.call(rbind, metrics)
rownames(metrics_FE) <- names(results_FE)
metrics_FE

metrics_OLS <- do.call(rbind, metrics_OLS)
rownames(metrics_OLS) <- names(results_OLS)
metrics_OLS


##Un-balance
###################***OLS########
#library(dplyr)
#library(plm)
#library(dplyr)
#
#set.seed(123)
#n <- 5000 
#alpha <- 1
#gamma <- -1
#beta1 <- 2
#beta2 <- -2
#theta <- 1
#delta <- 5 
#pR2 <- 0.7
#n_simulations=100
#model_results_OLS_random <- vector("list", n_simulations)
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
#      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
#      D = ifelse(Dtrue , 1, 0) 
#      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
#                           ifelse(time == 0 & D == 1, U_0_1,
#                                  ifelse(time == 1 & D == 0, U_1_0,
#                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
#      
#      # Simulate outcomes Y based on potential outcomes
#      Y <- alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected
#    })
#  
#  data_pR2_sim <- data  
#  D0 <- which(data_pR2_sim$D == 0)  
#  D1 <- which(data_pR2_sim$D == 1)  
#  
#  #Randomly remove some data to form a new non-response data set
#  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
#  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
#  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
#  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
#  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
#  data_pR2_sim <- data_pR2_sim[-rows_to_remove, ] 
#  
#  #  OLS regression
#  data_pR2_sim <- data_pR2_sim %>% mutate(did = d * time)
#  didOLS_random <- lm(Y ~ d + time + did, data = data_pR2_sim)
#  model_results_OLS_random[[i]] <- coef(didOLS_random)
#}
#
##----results---#
#results_matrix_OLS_random <- do.call(rbind, model_results_OLS_random)
#mean_coeffs_OLS_random <- apply(results_matrix_OLS_random, 2, mean)
#sd_coeffs_OLS_random <- apply(results_matrix_OLS_random , 2, sd)
#bias_OLS_random <- mean_coeffs_OLS_random - true_coeffs_OLS
#rmse_OLS_random <- sqrt(apply((results_matrix_OLS_random - true_coeffs_OLS)^2, 2, mean))
#
##----Print the results---#
#print(mean_coeffs_OLS_random)
#print(sd_coeffs_OLS_random)
#print(bias_OLS_random)
#print(rmse_OLS_random)
#
#
#
#
#
#############################################################
#
######***Random non-response-balance
#
############################################################
#
##balance
###################OLS########
#set.seed(123)
#model_results_OLS_random_bal <- vector("list", n_simulations)
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
#      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
#      D = ifelse(Dtrue , 1, 0) 
#      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
#                           ifelse(time == 0 & D == 1, U_0_1,
#                                  ifelse(time == 1 & D == 0, U_1_0,
#                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
#      
#      # Simulate outcomes Y based on potential outcomes
#      Y <- alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected
#    })
#  
#  data_pR2_sim <- data  
#  D0 <- which(data_pR2_sim$D == 0)  
#  D1 <- which(data_pR2_sim$D == 1)  
#  
#  #Randomly remove some data to form a new non-response data set
#  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
#  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
#  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
#  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
#  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
#  data_pR2_sim <- data_pR2_sim[-rows_to_remove, ] 
#  
#  # Create a balanced sample
#  ids_in_both_periods <- intersect(
#    data_pR2_sim %>% filter(time == 0) %>% select(id) %>% unlist(),
#    data_pR2_sim %>% filter(time == 1) %>% select(id) %>% unlist()
#  )
#  balanced_data <- data_pR2_sim %>% filter(id %in% ids_in_both_periods)
#  
#  # OLS regression on balanced data
#  balanced_data <- balanced_data %>% mutate(did = d * time)
#  didOLS_random <- lm(Y ~ d + time + did, data = balanced_data)
#  model_results_OLS_random_bal[[i]] <- coef(didOLS_random)
#}
#  
#
##----results---#
#results_matrix_OLS_random_bal <- do.call(rbind, model_results_OLS_random_bal)
#mean_coeffs_OLS_random_bal <- apply(results_matrix_OLS_random_bal, 2, mean)
#sd_coeffs_OLS_random_bal <- apply(results_matrix_OLS_random_bal , 2, sd)
#bias_OLS_random_bal <- mean_coeffs_OLS_random_bal - true_coeffs_OLS
#rmse_OLS_random_bal <- sqrt(apply((results_matrix_OLS_random_bal - true_coeffs_OLS)^2, 2, mean))
#
##----Print the results---#
#print(mean_coeffs_OLS_random_bal)
#print(sd_coeffs_OLS_random_bal)
#print(bias_OLS_random_bal)
#print(rmse_OLS_random_bal)
#
##balance
###################***FE########
#set.seed(123)
#model_results_FE_random_bal <- vector("list", n_simulations)
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
#      Dtrue  = time> 0 & gamma + beta1 + beta2 + theta + delta * C + E > 0  
#      D = ifelse(Dtrue , 1, 0) 
#      U_selected <- ifelse(time == 0 & D == 0, U_0_0,
#                           ifelse(time == 0 & D == 1, U_0_1,
#                                  ifelse(time == 1 & D == 0, U_1_0,
#                                         ifelse(time == 1 & D == 1, U_1_1, NA))))
#      
#      # Simulate outcomes Y based on potential outcomes
#      Y <- alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected
#    })
#  
#  data_pR2_sim <- data  
#  D0 <- which(data_pR2_sim$D == 0)  
#  D1 <- which(data_pR2_sim$D == 1)  
#  
#  #Randomly remove some data to form a new non-response data set
#  n_rows_to_remove_D0 <- ceiling(pR2 * length(D0))  
#  n_rows_to_remove_D1 <- ceiling(pR2 * length(D1))  
#  rows_to_remove_D0 <- sample(D0, size = n_rows_to_remove_D0, replace = FALSE)  
#  rows_to_remove_D1 <- sample(D1, size = n_rows_to_remove_D1, replace = FALSE)  
#  rows_to_remove <- unique(c(rows_to_remove_D0, rows_to_remove_D1))  
#  data_pR2_sim <- data_pR2_sim[-rows_to_remove, ] 
#  
#  # Create a balanced sample
#  ids_in_both_periods <- intersect(
#    data_pR2_sim %>% filter(time == 0) %>% select(id) %>% unlist(),
#    data_pR2_sim %>% filter(time == 1) %>% select(id) %>% unlist()
#  )
#  balanced_data <- data_pR2_sim %>% filter(id %in% ids_in_both_periods)
#  
#  # FE regression on balanced data
#  balanced_data <- balanced_data %>% mutate(did = d * time)
#  balanced_data_P <- pdata.frame(balanced_data, index=c("id", "time"))
#  didFE_random <- plm(Y ~ d + time + did, data = balanced_data_P, model="within")
#  model_results_FE_random_bal[[i]] <- coef(didFE_random)
#}
#
#
##----results---#
#results_matrix_FE_random_bal <- do.call(rbind, model_results_FE_random_bal)
#mean_coeffs_FE_random_bal <- apply(results_matrix_FE_random_bal, 2, mean)
#sd_coeffs_FE_random_bal <- apply(results_matrix_FE_random_bal , 2, sd)
#bias_FE_random_bal <- mean_coeffs_FE_random_bal - true_coeffs_FE
#rmse_FE_random_bal <- sqrt(apply((results_matrix_FE_random_bal - true_coeffs_FE)^2, 2, mean))
#
##----Print the results---#
#print(mean_coeffs_FE_random_bal)
#print(sd_coeffs_FE_random_bal)
#print(bias_FE_random_bal)
#print(rmse_FE_random_bal)
#
#
#
#