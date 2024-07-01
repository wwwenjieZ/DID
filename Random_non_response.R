

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

#data set simulation check

##############################################
#pR2=0.3
#pR2=0.4
pR2=0.15
n_rows_to_remove <- ceiling(pR2 * nrow(data))
rows_to_remove <- sample(seq_len(nrow(data)), size = n_rows_to_remove)
data_potnetial <- data[-rows_to_remove, ]
nrow(data_potnetial)

# Create a balanced sample
ids_in_both_periods <- intersect(
  data_potnetial %>% filter(time == 0) %>% select(id) %>% unlist(),
  data_potnetial %>% filter(time == 1) %>% select(id) %>% unlist()
)

data_potnetial_balance<-data_potnetial %>% filter(id %in% ids_in_both_periods)
nrow(data_potnetial_balance)

###############################################

#unbalance data set simulation 

##############################################
#pR2=0.3
#pR2=0.5
#pR2=0.15
pR2=0.05
n_simulations=1000
set.seed(123)
model_results_OLS_unbal <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2))
  
  data <- data %>%
    mutate(
      C = rep(rnorm(n, mean = 0, sd = 1), each = 2),
      U_0_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_0_1 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_1 = rnorm(n * 2, mean = 0, sd = 1),
      E = rnorm(n * 2, mean = 0, sd = 1)
    )%>%
    group_by(id) %>%
    mutate(
      Dtrue = any((time == 1) & (gamma + beta1 + beta2 + theta + delta * C + E > 0)),
      D = ifelse(Dtrue, 1, 0)
    ) %>%
    ungroup() %>%
    mutate(
      U_selected = case_when(
        time == 0 & d == 0 ~ U_0_0,
        time == 0 & d == 1 ~ U_0_1,
        time == 1 & d == 0 ~ U_1_0,
        time == 1 & d == 1 ~ U_1_1),
      Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected, did = d * time)
  
  n_rows_to_remove <- ceiling(pR2 * nrow(data))
 rows_to_remove <- sample(seq_len(nrow(data)), size = n_rows_to_remove)
  data<- data[-rows_to_remove, ]
  
  
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
      C = rep(rnorm(n, mean = 0, sd = 1), each = 2),
      U_0_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_0_1 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_1 = rnorm(n * 2, mean = 0, sd = 1),
      E = rnorm(n * 2, mean = 0, sd = 1)
    )%>%
    group_by(id) %>%
    mutate(
      Dtrue = any((time == 1) & (gamma + beta1 + beta2 + theta + delta * C + E > 0)),
      D = ifelse(Dtrue, 1, 0)
    ) %>%
    ungroup() %>%
    mutate(
      U_selected = case_when(
        time == 0 & d == 0 ~ U_0_0,
        time == 0 & d == 1 ~ U_0_1,
        time == 1 & d == 0 ~ U_1_0,
        time == 1 & d == 1 ~ U_1_1),
      Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected, did = d * time)

  
  #Randomly remove some data to form a new non-response data set
  n_rows_to_remove <- ceiling(pR2 * nrow(data))
  rows_to_remove <- sample(seq_len(nrow(data)), size = n_rows_to_remove)
  data<- data[-rows_to_remove, ]
  
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


# print result
metrics_FE <- do.call(rbind, metrics)
rownames(metrics_FE) <- names(results_FE)
metrics_FE

set.seed(123)
model_results_OLS <- vector("list", n_simulations)
for (i in 1:n_simulations) {
  data <- data.frame(
    id = rep(1:n, each = 2),   
    time = rep(0:1, times = n),  
    d = rep(sample(c(0, 1), n, replace = TRUE), each = 2))
  
  data <- data %>%
    mutate(
      C = rep(rnorm(n, mean = 0, sd = 1), each = 2),
      U_0_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_0_1 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_0 = rnorm(n * 2, mean = 0, sd = 1),
      U_1_1 = rnorm(n * 2, mean = 0, sd = 1),
      E = rnorm(n * 2, mean = 0, sd = 1)
    )%>%
    group_by(id) %>%
    mutate(
      Dtrue = any((time == 1) & (gamma + beta1 + beta2 + theta + delta * C + E > 0)),
      D = ifelse(Dtrue, 1, 0)
    ) %>%
    ungroup() %>%
    mutate(
      U_selected = case_when(
        time == 0 & d == 0 ~ U_0_0,
        time == 0 & d == 1 ~ U_0_1,
        time == 1 & d == 0 ~ U_1_0,
        time == 1 & d == 1 ~ U_1_1),
      Y = alpha + beta1 * time + beta2 * D + theta * (time * d) + C * D + U_selected, did = d * time)
  
  
  #Randomly remove some data to form a new non-response data set
  n_rows_to_remove <- ceiling(pR2 * nrow(data))
  rows_to_remove <- sample(seq_len(nrow(data)), size = n_rows_to_remove)
  data<- data[-rows_to_remove, ]

  
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



metrics_OLS <- do.call(rbind, metrics_OLS)
rownames(metrics_OLS) <- names(results_OLS)
metrics_OLS
metrics_FE

