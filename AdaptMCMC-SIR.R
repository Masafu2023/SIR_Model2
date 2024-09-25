# trial 2
# susceptible infected recovered (SIR)
library(deSolve)
library(tidyverse)
library(coda)
library(MCMCvis)

SIR_model <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  N <- S + I + R 
  
  beta <- parameters["beta"] # force of infection
  gamma <- parameters["gamma"] # recovery rate
  
  dS <- -beta * S * I / N 
  dI <- beta * S * I / N - gamma * I
  dR <- gamma * I
  
  list(c(dS, dI, dR))
}

# Parameters
beta <- 0.3      # Force of infection (FOI)
gamma <- 0.1     # Recovery rate == arbitrary values

# Initial conditions
initial_state <- c(S = 99, I = 1, R = 0)

# Time sequence (one year with daily time steps)
times <- seq(0, 365, by = 1)

#SIMULATION OF DATASET

# Run the SIR model
parameters <- c(beta = beta, gamma = gamma)
simulated_data <- ode(y = initial_state, times = times, func = SIR_model, parms = parameters) %>% data.frame()

# Adding a Time column for plotting
simulated_data$Time <- times
  
  

# Reshaping data to long format
simulated_data_long <- simulated_data %>%
  pivot_longer(cols = c(S, I, R), names_to = "State", values_to = "Count")

# Plotting all states on one graph
ggplot(simulated_data_long, aes(x = Time, y = Count, color = State)) +
  geom_line() +
  labs(title = "SIR Model Simulation", x = "Time (days)", y = "Count") +
  scale_color_manual(values = c("blue", "green", "red")) +
  theme_minimal()



# Log likelihood 
log_likelihood <- function(params, observed_data, initial_state, times) {
  # Parameters with the current estimate of beta (FOI)
  parameters <- c(params["beta"], gamma = 1/10) # Fixing gamma
  
  # Simulating the model
  out <- ode(y = initial_state, times = times, func = SIR_model, parms = parameters)
  model_data <- as.data.frame(out)
  
  # Calculating the log likelihood
  I_obs <- ceiling(observed_data$I)
  I_model <- model_data$I
  
  S_obs <- ceiling(observed_data$S)
  S_model <- model_data$S
  
  R_obs <- ceiling(observed_data$R)
  R_model <- model_data$R
  
  # Assuming Poisson distribution for number of infected individuals
  log_likelihood <- sum(dpois(I_obs, lambda = I_model, log = TRUE)) +
    sum(dpois(S_obs, lambda = S_model, log = TRUE))
  
  return(log_likelihood)
  
  # print(log_likelihood)
}


# MCMC settings
n_iter <- 5000
beta_init <-  c(0.1, 0.15, 0.2)
n_chains <- 3
beta_chain <- numeric(n_iter)
# beta_chains[1,] <- beta_init # initial values for each chain ...included to fit R-hat
sd_prop <- 0.001  # Initial standard deviation for proposal distribution
target_accept_rate <- 0.234  # Target acceptance rate for adaptive MCMC
adapt_rate <- 0.01  # Rate of adaptation for the proposal sd


# incorporating RHAT LOGIC
#storage for multiple chains
beta_chains <- matrix(NA, ncol = n_chains, nrow = n_iter)

#initializing chains with different initial values
beta_chains[1, ]<- beta_init

print(beta_chains[1, ])


# Prior distribution: Beta(2, 2)
prior <- function(beta) {
  return(dbeta(beta, 2, 2, log = TRUE))
}

loglik_curr <- log_likelihood(params=c(beta = beta_chain[1]), observed_data = simulated_data, 
                              initial_state = initial_state, times = times) + prior(beta_chain[1])

for (chain in 1:n_chains) {
  
  acceptance_counter <- 0 #Reset acceptance counter for each chain
  
  loglik_curr <- log_likelihood(params = c(beta = beta_chains[1, chain]),
                                
                                observed_data = simulated_data,
                                initial_state = initial_state,
                                times = times
                                )+
    prior(beta_chains[1,chain])
  
  
  
  for (i in 2:n_iter) {
    # Propose new beta from a normal distribution centered around the current value
    beta_proposed <- rnorm(1, mean = beta_chains[i - 1, chain], sd = sd_prop)
    
    if (beta_proposed > 0) {  # Ensure beta is positive
      loglik_prop <- log_likelihood(c(beta = beta_proposed), observed_data = simulated_data, 
                                    initial_state = initial_state, times = times) + 
        prior(beta_proposed)
    } else {
      loglik_prop <- -1E6
    }
    
    # Calculating acceptance probability
    acceptance_prob <- loglik_prop - loglik_curr
    
    # Metropolis-Hastings acceptance step
    if (log(runif(1)) < acceptance_prob) {
      beta_chains[i, chain] <- beta_proposed
      loglik_curr <- loglik_prop
      acceptance_counter <- acceptance_counter + 1
    } else {
      beta_chains[i, chain] <- beta_chains[i - 1, chain]
    }
    
    # Adaptive adjustment of proposal standard deviation (sd_prop)
    if (i > 100) {  # Adapt every 100 iterations
      acceptance_rate <- acceptance_counter / i
      sd_prop <- sd_prop * exp(adapt_rate * (acceptance_rate - target_accept_rate))
    }
  }
}




# Ploting the MCMC trace and histogram
burnin = 1000
mcmc_out <- as.mcmc(beta_chain[burnin:length(beta_chain)])

# Convert to data frame for ggplot
mcmc_df <- data.frame(Iteration = burnin:length(beta_chain), Beta = beta_chain[burnin:length(beta_chain)])



# Plot MCMC Trace for Beta
ggplot(mcmc_df, aes(x = Iteration, y = Beta)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "MCMC Trace for Beta (Adaptive MCMC)", x = "Iteration", y = "Beta") +
  theme_minimal()

parameters_fitted <- c(beta=0.276, gamma=0.1)
fitted_data <- ode(y = initial_state, times = times, func = SIR_model, parms = parameters_fitted) %>% data.frame()

# Convert simulated and fitted data to data frames
simulated_data_df <- as.data.frame(simulated_data)
fitted_data_df <- as.data.frame(fitted_data)

ggplot() +
  geom_line(data = simulated_data_df, aes(x = times, y = I), color = "blue", linewidth = 1) +
  geom_line(data = fitted_data_df, aes(x = times, y = I), color = "red", linewidth = 1, linetype = "dashed") +
  labs(title = "Simulated vs Fitted Infected Data", x = "Time (days)", y = "Infected Individuals") +
  theme_minimal()

# Histogram of posterior distribution using ggplot2
ggplot(data.frame(Beta = beta_chain), aes(x = Beta)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Posterior Distribution of Beta", x = "Beta", y = "Frequency") +
  theme_minimal()

beta_estimate <- mean(beta_chain)
cat("Estimated FOI (Beta):", beta_estimate, "\n")
cat("True FOI (Beta):", parameters["beta"], "\n")







# Rcpp code
library(Rcpp)
library(ggplot2)

# C++ SIR model function using Rcpp
cppFunction('
NumericMatrix SIR_model(double beta, double gamma, double S, double I, double R, int time_steps) {
    NumericMatrix results(time_steps, 3);
    
    for (int t = 0; t < time_steps; t++) {
        double N = S + I + R;
        double dS = -beta * S * I / N;
        double dI = beta * S * I / N - gamma * I;
        double dR = gamma * I;
        
        S += dS;
        I += dI;
        R += dR;
        
        results(t, 0) = S;
        results(t, 1) = I;
        results(t, 2) = R;
    }
    
    return results;
}
')

# Parameters and initial conditions for Rcpp model
beta <- 0.3
gamma <- 0.1
initial_state <- c(S = 99, I = 1, R = 0)
time_steps <- 366

# Simulating the data using Rcpp SIR model
simulated_data <- SIR_model(beta, gamma, initial_state[1], initial_state[2], initial_state[3], time_steps)

# Converting the simulated data
simulated_data_df <- as.data.frame(simulated_data)
colnames(simulated_data_df) <- c("S", "I", "R")
simulated_data_df$time <- 1:time_steps

# Reshaping the data to long format for ggplot2 (so S, I, R are different lines)
simulated_data_long <- simulated_data_df %>%
  pivot_longer(cols = c("S", "I", "R"), names_to = "Compartment", values_to = "Count")

# Plots; the three lines (S, I, R) on the same graph
ggplot(simulated_data_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line(size = 1) +
  labs(title = "SIR Model: Susceptible, Infected, Recovered (Rcpp)", x = "Time (days)", y = "Population Count") +
  theme_minimal()






# MCMC-VISUALIZATION using MCMCvis;
library(MCMCvis)
mcmc_matrix <- as.matrix(mcmc_df)
chain1 <- as.mcmc(mcmc_matrix[1:1500, ])  # First half of the rows
chain2 <- as.mcmc(mcmc_matrix[1501:3000, ])  # Second half of the rows
mcmc_list <- mcmc.list(chain1, chain2)
MCMCsummary(mcmc_list)
MCMCsummary(mcmc_list, 
            params = 'all', 
            ISB = FALSE, 
            exact = FALSE, 
            round = 2)

#MCMCpstr
MCMCpstr(mcmc_list, 
         params = 'Beta',
         func = mean,
         # round = 2,
         type = 'summary')

ex <- MCMCpstr(mcmc_list, type = 'chains')
dim(ex$Beta)


# MCMCtrace
MCMCtrace(mcmc_list, 
          params = 'Beta', 
          type = 'density',
          ISB = FALSE, 
          exact = TRUE,
          ind = TRUE,
          pdf = FALSE,
          iter = 100, 
          open_pdf = FALSE,
          # priors = beta,
          filename = 'Mypdf',
          plot = FALSE,
          PPO_out = TRUE,
          post_zm = FALSE,
          wd = "DIRECTORY_HERE")

# MCMCchains
ex <- MCMCchains(mcmc_list, params = 'Beta')




####################

# Variable to track acceptance rate
acceptance_counter <- 0

# Run the adaptive MCMC
for (i in 2:n_iter) {
  # Propose new beta from a normal distribution centered around the current value
  beta_proposed <- rnorm(1, mean = beta_chain[i-1], sd = sd_prop)
  
  if (beta_proposed > 0) {  # Ensure beta is positive
    loglik_prop <- log_likelihood(c(beta = beta_proposed), observed_data = simulated_data, 
                                  initial_state = initial_state, times = times) + prior(beta_proposed)
  } else {
    loglik_prop <- -1E6
  }
  
  # Calculating acceptance probability
  acceptance_prob <- loglik_prop - loglik_curr
  
  # Metropolis-Hastings acceptance step
  if (log(runif(1)) < acceptance_prob) {
    beta_chain[i] <- beta_proposed
    loglik_curr <- loglik_prop
    acceptance_counter <- acceptance_counter + 1
  } else {
    beta_chain[i] <- beta_chain[i-1]
  }
  
  # Adaptive adjustment of proposal standard deviation (sd_prop)
  if (i >100) {  # Adapt every 100 iterations
    acceptance_rate <- acceptance_counter / i
    # Update proposal sd based on acceptance rate
    sd_prop <- sd_prop * exp(adapt_rate * (acceptance_rate - target_accept_rate))
    # acceptance_counter <- 0  # Reset counter
  }
}