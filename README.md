In each R file, I repeated the following steps

First, we creates a data frame with columns id, time, and adding variables and simulating outcomes:
  C, U_0_0, U_0_1, U_1_0, U_1_1, E: Random normal variables representing different components of the model.
  D_true: Defines the true treatment indicator based on time.(time > 0 & (gamma + beta1 + beta2 + theta + delta * C + E) > 0)
  U_selected: Chooses the appropriate potential outcome based on D and time.
  Y: Simulates the outcome variable using a linear combination of the defined variables and parameters.

Second, we adds the post and did interaction terms and estimates the model using OLS (lm) and fixed effects (plm) methods directly for extracting the coefficients for further comparison.

Third, Simulates non-response(**4 senerios**)
     (1) randon non-response
     (2) non-response in post treatment group 
     (3) non-response in pre treatment group 
     (4) asymmetric selective non-response in both of group 
     In senerio (1), randon non-response probabiity is pR2= 0.2, for senerios(2) and (3), we Simulating Non-Response Probability by calculates the cumulative probability (CDF) of the standard normal distribution for each value of Y, if the non_response_probability is greater than 0.5, indicating non-response.( With the probability set to 0.5, the missing data in the unbalnce data set is around two thousand, which is roughly the same as scenario (1), because we don't want the difference in results between them to be due to how much is missing)  In higher non-response example, pR2= 0.5 , non_response_probability is o.1 for senerios (2) and (3), and 0.5 for senerio (4).

Fourth, simulation for Unbalanced Data Set, we repeats the data generation and OLS estimation 1000 times to account for variability in unbalanced data, and alculating bias, standard deviation, and RMSE for the OLS estimates. 

Fifth, Repeats the data generation and FE/OLS estimation 1000 times for balanced data. Calculates bias, standard deviation, and RMSE for the FE/OLS estimates.
