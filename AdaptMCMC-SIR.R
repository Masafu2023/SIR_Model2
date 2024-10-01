---
title: "ConvergentMCMC_SIR"
author: "Brian Masafu"
date: "2024-09-26"
output:
  pdf_document: default
  html_document: default
---


# SIR_CompartmentalWith_3_AgeStructuredComparments 
# libraries
```{r}
library(deSolve)
library(tidyverse)
library(coda)
library(MCMCvis)
library(lattice)

```


# SIR model function with age compartments
```{r}
SIR_model_age <- function(time, state, parameters) {
  # Children (age group 1; 0-17 years)
  S1 <- state[1]
  I1 <- state[2]
  R1 <- state[3]
  
  # Adults (age group 2; 18-59 years)
  S2 <- state[4]
  I2 <- state[5]
  R2 <- state[6]
  
  # Elderly (age group 3; 60 - 99 years)
  S3 <- state[7]
  I3 <- state[8]
  R3 <- state[9]
  
  # Total population for each age group
  N1 <- S1 + I1 + R1
  N2 <- S2 + I2 + R2
  N3 <- S3 + I3 + R3
  
  # Transmission rates (force of infection) for each age group & recovery rate
  beta1 <- parameters["beta1"]
  beta2 <- parameters["beta2"]
  beta3 <- parameters["beta3"]
  gamma <- parameters["gamma"]
  
  # Birth and death rates
  birth_rate <- parameters["birth_rate"]
  death_rate <- parameters["death_rate"]
  
  # Vaccine coverage and efficacy (all-or-nothing or leaky immunity)
  vax_coverage <- parameters["vax_coverage"]
  vax_efficacy <- parameters["vax_efficacy"]
  
  # Aging rate (children becoming adults, adults becoming elderly)
  aging_rate1to2 <- parameters["aging_rate1to2"]
  aging_rate2to3 <- parameters["aging_rate2to3"]
  
  # Equations for age group 1 (Children)
  dS1 <- -beta1 * S1 * I1 / N1 + birth_rate * (1 - vax_coverage) * N1 - death_rate * S1 - aging_rate1to2 * S1
  dI1 <- beta1 * S1 * I1 / N1 - gamma * I1 - death_rate * I1 - aging_rate1to2 * I1
  dR1 <- gamma * I1 + birth_rate * vax_coverage * vax_efficacy * N1 - death_rate * R1 - aging_rate1to2 * R1
  
  # Equations for age group 2 (Adults)
  dS2 <- -beta2 * S2 * I2 / N2 + aging_rate1to2 * S1 - death_rate * S2 - aging_rate2to3 * S2
  dI2 <- beta2 * S2 * I2 / N2 + aging_rate1to2 * I1 - gamma * I2 - death_rate * I2 - aging_rate2to3 * I2
  dR2 <- gamma * I2 + aging_rate1to2 * R1 - death_rate * R2 - aging_rate2to3 * R2
  
  # Equations for age group 3 (Elderly)
  dS3 <- -beta3 * S3 * I3 / N3 + aging_rate2to3 * S2 - death_rate * S3
  dI3 <- beta3 * S3 * I3 / N3 + aging_rate2to3 * I2 - gamma * I3 - death_rate * I3
  dR3 <- gamma * I3 + aging_rate2to3 * R2 - death_rate * R3
  
  # Returning the rates of change for all compartments
  list(c(dS1, dI1, dR1, dS2, dI2, dR2, dS3, dI3, dR3))
}

```


# Parameters
```{r}
# Define parameters
beta1 <- 0.4   # Transmission rate for children
beta2 <- 0.3   # Transmission rate for adults
beta3 <- 0.25  # Transmission rate for elderly
gamma <- 0.1   # Recovery rate (same for all age groups)
birth_rate <- 0.01  # Birth rate for youngest group
death_rate <- 0.01  # Death rate 
vax_coverage <- 0.6  # 60% vaccination coverage at birth
vax_efficacy <- 0.8  # 80% vaccine efficacy (all-or-nothing)
aging_rate1to2 <- 0.02  # Aging rate from children to adults
aging_rate2to3 <- 0.015 # Aging rate from adults to elderly

# Parameters vector
parameters <- c(beta1 = beta1, beta2 = beta2, beta3 = beta3, gamma = gamma, 
                birth_rate = birth_rate, death_rate = death_rate, 
                vax_coverage = vax_coverage, vax_efficacy = vax_efficacy,
                aging_rate1to2 = aging_rate1to2, aging_rate2to3 = aging_rate2to3)

```


# Initial conditions for each age group &time sequence
```{r}
# Initial conditions for each age group
initial_state <- c(S1 = 99, I1 = 1, R1 = 0,  # Children
                   S2 = 99, I2 = 1, R2 = 0,  # Adults
                   S3 = 99, I3 = 1, R3 = 0)  # Elderly

# Time sequence (one year with daily time steps)
times <- seq(0, 365, by = 1)

```

# Model simulation and plots
```{r}
# Simulating the SIR model with age structure
simulated_data <- ode(y = initial_state, times = times, func = SIR_model_age, parms = parameters) %>% data.frame()

# observed data columns as numeric
simulated_data[, 2:10] <- lapply(simulated_data[, 2:10], as.numeric)

#  Time column for plotting
simulated_data$Time <- times

# Reshaping data to long format for plotting
simulated_data_long <- simulated_data %>%
  pivot_longer(cols = -Time, names_to = "Compartment", values_to = "Count")

# Ploting the SIR simulation 
ggplot(simulated_data_long, aes(x = Time, y = Count, color = Compartment)) +
  geom_line() +
  labs(title = "SIR Model Simulation with Birth, Aging, Death, and Vaccination", 
       x = "Time (days)", y = "Count") +
  scale_color_manual(values = c("blue", "green", "red",      # Children
                                "blue4", "green4", "red4",  # Adults
                                "blue2", "green2", "red2",  # Elderly
                                "purple")) +  # Additional color to make it 10
  coord_cartesian(ylim = c(0, 100)) +  # Set y-axis limits from 0 to 100
  theme_minimal()




```
SIR_plots showing different levels of force of infections and recovery rates among the 3 age structured population with (Susceptible, Infected and Recovered States)

# Log-likelihood function 
-Measures how likely is it to observe the given data for different values of the models parameters.
($\text{P}(\text{data}|\theta)$)
$\theta$ = set of parameters for the model
data = observed data.
In Bayesian statistics, likelihood combines with the prior distribution to update the parameters resulting in posterior distribution.


```{r}
log_likelihood_age <- function(params, observed_data, initial_state, times) {
  parameters <- c(params["beta1"], params["beta2"], params["beta3"], params["death_rate"], 
                  params["birth_rate"], params["vax_efficacy"], params["vax_coverage"], gamma = 1/10) 
  
  out <- ode(y = initial_state, times = times, func = SIR_model_age, parms = parameters)
  
  if (any(params <= 0)) {
        return(-Inf)  # Reject invalid parameters
  }
  
  model_data <- as.data.frame(out)
  observed_data[, 2:10] <- lapply(observed_data[, 2:10], as.numeric)
  
  log_likelihood <- 0
  for (col in 2:10) {
    obs <- ceiling(observed_data[[col]])
    model <- model_data[[col]]
    
    if (any(is.na(obs)) || any(is.na(model))) {
      return(-Inf)
    }
    
    log_likelihood <- log_likelihood + sum(dpois(obs, lambda = model, log = TRUE))
  }
  
  return(log_likelihood)
}





```


# MCMC settings
- used to estimate posterior distribution
Markov chain - is a sequence of random variables where each variable(or state) depends on the previous one 
Monte Carlo - Estimates distribution of interest.
($\text{P}(X_{t}+ _1| \text{X}_{t1}, X_{t2}....X_{t}_{t} = \text{P}(X_{t}+_1|X_{t}$)
```{r}
n_iter <- 5000
beta_init <- c(0.1, 0.15, 0.2)  # Initial beta values for MCMC
n_chains <- 3 # Number of chains
sd_prop <- rep(0.001, n_chains)  # Initial SD for each chain
target_accept_rate <- 0.234  # Target acceptance rate for adaptive MCMC
adapt_rate <- 0.01  # Rate of adaptation for the proposal SD


```

# Incorporating RHAT logic: storage for multiple chains& initializing the chain
Rhat(potential scale reduction factor) - Measures whether or not an MCMC algorithm converged. It checks the distribution of a chain (after warm up ) the same as the distribution of the second half of the chain.
2. If the algorithm starts at 2 different places and chain left to warm up, both the chains have the same distribution.

```{r}
# Incorporating RHAT logic: storage for multiple chains
beta_chains <- matrix(NA, ncol = n_chains, nrow = n_iter)

# Initializing chains with different initial values
beta_chains[1, ] <- beta_init

# Prior distribution: Beta(2, 2)
prior <- function(beta) {
    if (beta <= 0 || beta >= 1) {
        return(-Inf)  # Reject invalid beta
    }
    return(dbeta(beta, 2, 2, log = TRUE))  # Return log density
}
```

# MCMC loop for each chain
```{r}

for (chain in 1:n_chains) {
  
  acceptance_counter <- 0  # Resetting acceptance counter for each chain
  loglik_curr <- log_likelihood_age(params = c(beta1 = beta_chains[1, chain], 
                                               beta2 = beta_chains[1, chain], 
                                               beta3 = beta_chains[1, chain], 
                                               birth_rate = beta_chains[1, chain], 
                                               death_rate = beta_chains[1, chain], 
                                               vax_coverage = beta_chains[1, chain], 
                                               vax_efficacy = beta_chains[1, chain]),
                                    observed_data = simulated_data,
                                    initial_state = initial_state,
                                    times = times) + prior(beta_chains[1, chain])
  
  if (is.nan(loglik_curr)) {
    loglik_curr <- -Inf
    print(paste("Warning: Initial log likelihood is NaN for chain", chain))
  }
  
  for (i in 2:n_iter) {
    beta_proposed <- rnorm(1, mean = beta_chains[i - 1, chain], sd = sd_prop[chain])
    
    # Check for NA or invalid proposed beta
    if (is.na(beta_proposed) || beta_proposed <= 0) {
      loglik_prop <- -Inf  # Penalize invalid proposals
      print(paste("Warning: Proposed beta is NA or non-positive at iteration", i, "in chain", chain))
    } else {
      loglik_prop <- log_likelihood_age(params = c(beta1 = beta_proposed, 
                                                   beta2 = beta_proposed, 
                                                   beta3 = beta_proposed,
                                                   birth_rate = beta_proposed, 
                                                   death_rate = beta_proposed, 
                                                   vax_coverage = beta_proposed, 
                                                   vax_efficacy = beta_proposed),
                                        observed_data = simulated_data,
                                        initial_state = initial_state,
                                        times = times) + prior(beta_proposed)
      
      if (is.nan(loglik_prop)) {
        loglik_prop <- -Inf
        print(paste("Setting loglik_prop to -Inf due to NaN at iteration", i, "in chain", chain))
      }
    }
    
    # Check loglik_curr before calculating acceptance probability
    if (is.nan(loglik_curr)) {
      loglik_curr <- -Inf
      print(paste("Setting loglik_curr to -Inf due to NaN at iteration", i, "in chain", chain))
    }
    
    # Compute acceptance probability only if valid log likelihood values
    if (!is.nan(loglik_curr) && !is.nan(loglik_prop)) {
      acceptance_prob <- loglik_prop - loglik_curr
      
      if (is.nan(acceptance_prob)) {
        print(paste("Warning: acceptance_prob is NaN at iteration", i, "in chain", chain))
      } else {
        # Metropolis-Hastings acceptance step
        if (log(runif(1, min = 1e-10, max = 1)) < acceptance_prob) {
          beta_chains[i, chain] <- beta_proposed
          loglik_curr <- loglik_prop
          acceptance_counter <- acceptance_counter + 1
        } else {
          beta_chains[i, chain] <- beta_chains[i - 1, chain]
        }
      }
    }
    
    # Adaptive adjustment of proposal standard deviation (sd_prop)
    if (i > 100) {  # Adapt every 100 iterations
      acceptance_rate <- acceptance_counter / i
      sd_prop[chain] <- sd_prop[chain] * exp(adapt_rate * (acceptance_rate - target_accept_rate))
      
      # Ensure sd_prop doesn't become negative or too small
      if (sd_prop[chain] < 1e-6) {
        sd_prop[chain] <- 1e-6  # Minimum bound to avoid degenerate proposals
      }
    }
  }
}


```


#converting beta chains

```{r}
burnin <- 1000
mcmc_out <- mcmc.list(
  as.mcmc(beta_chains[burnin:n_iter, 1, drop = FALSE]),
  as.mcmc(beta_chains[burnin:n_iter, 2, drop = FALSE]),
  as.mcmc(beta_chains[burnin:n_iter, 3, drop = FALSE])
)

# Converting beta_chains matrix into mcmc objects 
mcmc_chain1 <- as.mcmc(beta_chains[, 1])
mcmc_chain2 <- as.mcmc(beta_chains[, 2])
mcmc_chain3 <- as.mcmc(beta_chains[, 3])

# Combining all chains into an mcmc.list object
mcmc_combined <- mcmc.list(mcmc_chain1, mcmc_chain2, mcmc_chain3)

# Calculating R-hat using rs tan's summary function
rhat_values <- rstan::Rhat(as.matrix(mcmc_combined))

```


#printing R-hat values
```{r}
# Printing R-hat values for convergence diagnostics
print(rhat_values)
```
Rhat values close to 1(1.1) indicates that the chain has converged. Values greater than 1.1 suggests that the chain has not yet converged and you may need to run more iteration. `Rhat calculations` You need multiple chains to calculate Rhat.



# Geweke statistics
```{r}
library(coda)
# Applying  geweke.diag to each chain individually
geweke_results <- list(
  geweke.diag(mcmc_chain1),
  geweke.diag(mcmc_chain2),
  geweke.diag(mcmc_chain3)
)

#  structure of the results
str(geweke_results)

#  structure and class of geweke_results
str(geweke_results)

#  contents of geweke_results 
print(geweke_results)

# checking  its names if its a list
if (is.list(geweke_results)) {
  print(names(geweke_results))  # This will show the names of the list components
  # Check the structure of the first element
  str(geweke_results[[1]])
}


# Assuming geweke_results is a list
if (is.list(geweke_results)) {
  # Initialize empty vectors to hold parameters and z-values
  parameter_names <- c()
  z_values <- c()
  
  for (i in seq_along(geweke_results)) {
    # Checking if each component has the necessary information
    if ("z" %in% names(geweke_results[[i]])) {
      z_values <- c(z_values, geweke_results[[i]]$z)
      parameter_names <- c(parameter_names, rownames(geweke_results[[i]]))
    }
  }
  

  if (length(parameter_names) > 0 && length(z_values) > 0) {
    geweke_values <- data.frame(
      Parameter = parameter_names,
      z = z_values
    )
  } else {
    warning("No Z-values or parameter names found.")
  }
}
# If geweke_results is a data frame or matrix
if (is.data.frame(geweke_results) || is.matrix(geweke_results)) {
  parameter_names <- rownames(geweke_results)  # Assuming rownames hold parameter names
  z_values <- geweke_results[, "z"]  # Adjust based on actual column name
  
  # Creating the data frame
  geweke_values <- data.frame(
    Parameter = parameter_names,
    z = z_values
  )
}

```
based on z - value(-1.132, -0.9043, -1.027) the chain has converged. since z values close to 0 indicates that the chain has converged.    
This chains show little to no significant drift between the early and late samples, indicating good mixing and convergence of the MCMC process.
Availability of windows fractions despite no Z - values ,suggest that the model has likely converged.

#plots for trace mcmc
```{r}

#  trace plot using coda's plot function
plot(mcmc_combined, trace = TRUE, density = FALSE, main = "Trace Plots for MCMC Chains")

# Assigning column names to the beta_chains matrix 
colnames(beta_chains) <- c("Beta1", "Beta2", "Beta3")

# Converting each chain to mcmc objects again
mcmc_chain1 <- as.mcmc(matrix(beta_chains[, "Beta1"], ncol = 1, dimnames = list(NULL, "Beta")))
mcmc_chain2 <- as.mcmc(matrix(beta_chains[, "Beta2"], ncol = 1, dimnames = list(NULL, "Beta")))
mcmc_chain3 <- as.mcmc(matrix(beta_chains[, "Beta3"], ncol = 1, dimnames = list(NULL, "Beta")))

# Combining  into an mcmc.list object
mcmc_combined <- mcmc.list(mcmc_chain1, mcmc_chain2, mcmc_chain3)

# MCMCtrace 
MCMCtrace(mcmc_combined, params = "Beta", pdf = FALSE)

# Ploting posterior density using coda's densityplot function
densityplot(mcmc_combined, main = "Posterior Density of Beta")


# posterior  visualization using MCMCvis
MCMCplot(mcmc_combined, params = "Beta", main = "Posterior Distribution of Beta")

```

# Calculating Gelman and Rubin's diagnostic (R-hat)
```{r}
gelman.diag(mcmc_combined)
```
This MCMC chains have likely converged to the same target distribution, and there is little difference between the chains' behavior.Values close to 1 (typically less than 1.1) are considered evidence of convergence.

# MCMCvis to summarize the chains including R-hat and effective sample size
```{r}

MCMCsummary(mcmc_combined, Rhat = TRUE, n.eff = TRUE)
```

# rstan for comprehensive summary including n_eff and other statistics
```{r}

rstan::summary(mcmc_combined)
```

# MCMC Trace for Beta && # Histogram of posterior distribution
```{r}
# MCMC Trace for Beta
mcmc_df <- data.frame(Iteration = burnin:n_iter, Beta = c(beta_chains[burnin:n_iter, ]))
ggplot(mcmc_df, aes(x = Iteration, y = Beta)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "MCMC Trace for Beta (Adaptive MCMC)", x = "Iteration", y = "Beta") +
  theme_minimal()

# Histogram of posterior distribution
ggplot(data.frame(Beta = c(beta_chains[burnin:n_iter, ])), aes(x = Beta)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Posterior Distribution of Beta", x = "Beta", y = "Frequency") +
  theme_minimal()

```
