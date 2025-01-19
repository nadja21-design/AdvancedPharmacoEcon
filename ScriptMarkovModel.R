# Load packages
library(heemod)
library(dplyr)

# General death probabilities for UHR
death_prob = 0.003 # defined outside define_parameters() because there was a bug
rr = 0.367 # effectiveness of CBT compared to TAU from RCT 


# Step 1: Define parameters
param <- define_parameters(
  
  age_init = 14, # assuming our cohort starts at age 14
  
  # age increases with cycles
  age = age_init + model_time,
  
  # Transition probabilities 
  p_uh_to_ps = 0.2152,                          # UHR → Psychosis
  p_uh_to_ns = 0.4465,                          # UHR → No Symptoms
  p_uh_to_ps_cbt = p_uh_to_ps*rr,               # UHR → Psychosis for CBT group 
  p_ps_to_d = death_prob * 2.58 + 0.0057,       # Psychosis → Death
  p_ps_to_pp = 0.7862,                          # (markov paper) Psychosis → Post-Psychosis + remission
  #p_ps_to_pp = 1-p_ps_to_d, # no longer correct                     
  p_pp_to_ps = 1 - exp(-(-log(1 - 0.51) / 3)),  # Recurrence over 3 years converted to annual probability
  p_pp_to_d = death_prob + 0.0057,              # Death for Post-Psychosis 
  p_pp_to_pp = 1 - (p_pp_to_ps + 0.0057 + p_pp_to_d), # Complement
  p_ns_to_d = death_prob, 
  p_ns_to_ns = 1 - p_ns_to_d,
  p_uh_to_uh = 1 - (p_uh_to_ps + p_uh_to_ns + death_prob),  # Complement for staying in UHR, MIGHT NEED TO REMOVE
   
  # Discount rates
  dr_costs = 0.04,       # Costs discount rate
  dr_health = 0.015,     # Health effects discount rate
  
  # Costs
  cost_uh = 4874,         # Cost per cycle in UHR
  cost_ps = 7200,         # Cost per cycle in Psychosis
  cost_pp = 5019,         # Cost per cycle in Post-Psychosis
  cost_ns = 3970,         # Cost per cycle in No Symptoms, reminder to double check this parameter
  
  #Cost TAU
  cost_tau = 0, #NOTE TO CHANGE THIS ONE 
  
  #Cost CBT
  #cost_cbt = ifelse(model_time <= 1, 1924, 1000),
  cost_cbt = ifelse(model_time <= 1, 1924, 0),
  
  # Utilities
  util_uh = 0.64,        # Utility for UHR
  util_ps = 0.34,        # Utility for Psychosis
  util_pp = 0.362,       # Utility for Post-Psychosis
  util_ns = 0.8,         # Utility for No Symptoms
)

# Transition matrix for CBT
tm_cbt <- define_transition(
  C, p_uh_to_ps_cbt, 0, p_uh_to_ns, death_prob,   # UHR row
  0, 0, p_ps_to_pp, 0, C,                         # Psychosis row
  0, p_pp_to_ps, C, 0, p_pp_to_d,                 # Post-Psychosis row
  0, 0, 0, C, death_prob,                         # No Symptoms row
  0, 0, 0, 0, 1                                   # Death row
)

# Transition matrix for TAU
tm_tau <- define_transition(
  C, p_uh_to_ps, 0, p_uh_to_ns, death_prob,       # UHR row
  0, 0, p_ps_to_pp, 0, C,                         # Psychosis row
  0, p_pp_to_ps, C, 0, p_pp_to_d,                 # Post-Psychosis row
  0, 0, 0, C, death_prob,                         # No Symptoms row
  0, 0, 0, 0, 1                                   # Death row
)



# Step 3: Define health states CBT
cbt_state_uhr <- define_state(
  cost = discount(cost_uh + cost_cbt, dr_costs),
  utility = discount(util_uh, dr_health) 
)

cbt_state_ps <- define_state(
  cost = discount(cost_ps, dr_costs),
  utility = discount(util_ps, dr_health)
)
cbt_state_pp <- define_state(
  cost = discount(cost_pp, dr_costs),
  utility = discount(util_pp, dr_health)
)
cbt_state_ns <- define_state(
  cost = discount(cost_ns, dr_costs),
  utility = discount(util_ns, dr_health)
)
cbt_state_d <- define_state(
  cost = 0,
  utility = 0
)

# CBT strategy
strategy_cbt <- define_strategy(
  transition = tm_cbt,  # Use the defined transition matrix
  cbt_state_uhr,
  cbt_state_ps,
  cbt_state_pp,
  cbt_state_ns,
  cbt_state_d
)

# Step 3: Define health states TAU
tau_state_uhr <- define_state(
  cost = discount(cost_uh + cost_tau, dr_costs),
  utility = discount(util_uh, dr_health) 
)
tau_state_ps <- define_state(
  cost = discount(cost_ps, dr_costs),
  utility = discount(util_ps, dr_health)
)
tau_state_pp <- define_state(
  cost = discount(cost_pp, dr_costs),
  utility = discount(util_pp, dr_health)
)
tau_state_ns <- define_state(
  cost = discount(cost_ns, dr_costs),
  utility = discount(util_ns, dr_health)
)
tau_state_d <- define_state(
  cost = 0,
  utility = 0
)

# TAU strategy
strategy_tau <- define_strategy(
  transition = tm_tau,  # Use the defined transition matrix
  tau_state_uhr,
  tau_state_ps,
  tau_state_pp,
  tau_state_ns,
  tau_state_d
)

# Step 5: Run the Markov models
model_cbt <- run_model(
  strategy = strategy_cbt,
  parameters = param,
  cycles = 10,               # Number of cycles
  init = c(1000, 0, 0, 0, 0),
  cost = cost,
  effect = utility,
  method = "life-table"      # Use life-table method
)

model_tau <- run_model(
  strategy = strategy_tau,
  parameters = param,
  cycles = 10,               # Number of cycles
  init = c(1000, 0, 0, 0, 0),
  cost = cost,
  effect = utility, 
  method = "life-table"      # Use life-table method
)

# Summarize and visualize the results
summary(model_cbt)
plot(model_cbt)

summary(model_tau)
plot(model_tau)

# Combined model
results <- run_model(
  cbt = strategy_cbt,
  tau = strategy_tau,
  parameters = param,
  cycles = 10,
  init = c(1000, 0, 0, 0, 0),
  cost = cost,
  effect = utility, 
  method = "life-table"
)

summary(results)

plot(results)

# Adding probabilistic component 

psax <- define_psa(
  # Probabilistic costs
  cost_uh ~ gamma(mean = 4874, sd = sqrt(4874)),         # Cost per cycle in UHR
  cost_ps ~ gamma(mean = 7200, sd = sqrt(7200)),         # Cost per cycle in Psychosis
  cost_pp ~ gamma(mean = 5019, sd = sqrt(5019)),         # Cost per cycle in Post-Psychosis
  cost_ns ~ gamma(mean = 3970, sd = sqrt(3970)),         # Cost per cycle in No Symptoms
  #cost_tau ~ gamma(mean = 1000, sd = sqrt(1000)),        # Cost for TAU
  #cost_cbt ~ gamma(mean = 1924, sd = sqrt(1924)),        # Cost for CBT

  
  # Probabilistic utilities using beta distribution
  util_uh ~ beta(shape1 = 64, shape2 = 36),              # Utility for UHR (mean ~0.64)
  util_ps ~ beta(shape1 = 34, shape2 = 66),              # Utility for Psychosis (mean ~0.34)
  util_pp ~ beta(shape1 = 362, shape2 = 638),            # Utility for Post-Psychosis (mean ~0.362)
  util_ns ~ beta(shape1 = 80, shape2 = 20),              # Utility for No Symptoms (mean ~0.8)                    # No Symptoms

  # needs debugging
  # Probabilistic transitions using multinomial (with explicit naming)
  p_uh_to_uh + p_uh_to_ps + p_uh_to_ns + death_prob ~ multinomial(215, 447, 335, 3),   # UHR transitions
  p_ps_to_d ~ binomial(prob = 0.006, size = 1),                                # Psychosis transitions
  p_pp_to_pp + p_pp_to_ps + p_pp_to_d ~ multinomial(510, 362, 5),         # Post-Psychosis transitions
  p_ns_to_d ~ binomial(prob = 0.001, size = 1000)                         # No Symptoms transitions
  
  )

  

pm <- run_psa(
  model = results,
  psa = psax,
  N = 1000
)



summary(
  pm, 
  threshold = c(1000, 5000, 6000, 1e4))

### Result interpretation

# The results of the analysis can be plotted on the cost-effectiveness plane. 
# We can see there seem to be little uncertainty on the costs compared to the
# uncertainty on the effects, resulting in an uncertainty cloud that looks like a line.

plot(pm, type = "ce")

# And as cost-effectiveness acceptability curves or EVPI:

plot(pm, type = "ac", max_wtp = 10000, log_scale = FALSE)
plot(pm, type = "evpi", max_wtp = 10000, log_scale = FALSE)

# A covariance analysis can be performed on strategy results:

# plot(pm, type = "cov") doesn't work but might not be needed

# Or on the difference between strategies:

plot(pm, type = "cov", diff = TRUE, threshold = 5000)

# As usual plots can be modified with the standard ggplot2 syntax.

library(ggplot2)

plot(pm, type = "ce") +
  xlab("Life-years gained") +
  ylab("Additional cost") +
  scale_color_brewer(
    name = "Strategy",
    palette = "Set1"
  ) +
  theme_minimal()

#####

# Prepare PSA results for BCEA
# Extract costs and utilities for CBT and TAU
costs <- c(39318.47, 36365.76)
effects <- c(3.085737, 2.710254)
wtp <- seq(1000, 20000, by = 1000)  # Range of WTP thresholds
# NMB for each WTP threshold
nmb <- sapply(wtp, function(w) {
  cbt_nmb <- effects[1] * w - costs[1]  # NMB for CBT
  tau_nmb <- effects[2] * w - costs[2]  # NMB for TAU
  return(c(cbt_nmb, tau_nmb))
})

# Transpose the result for easier interpretation
nmb <- t(nmb)
colnames(nmb) <- c("CBT", "TAU")
# Calculate the EVPPI as the difference in NMB
evppi <- apply(nmb, 1, function(nmb_row) {
  max(nmb_row) - mean(nmb_row)  # Maximum NMB - average NMB
})

# Example: Define PSA results for parameters manually
psa_results <- list(
  costs = matrix(c(39318.47, 36365.76), nrow = 1, dimnames = list(NULL, c("cbt", "tau"))),
  effects = matrix(c(3.085737, 2.710254), nrow = 1, dimnames = list(NULL, c("cbt", "tau"))),
  p_uh_to_ps = c(0.2152, 0.22, 0.21),  # Example: Samples of UHR → Psychosis probabilities
  p_ps_to_pp = c(0.7862, 0.78, 0.79),  # Example: Samples of Psychosis → Post-Psychosis probabilities
  rr = c(0.367, 0.35, 0.38)            # Example: Samples of CBT relative risk
)

