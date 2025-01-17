# Load packages
library(heemod)
library(dplyr)

# General death probabilities for UHR
death_prob = 0.003 # defined outside define_parameters() because there was a bug
rr = 0.367 # effectiveness of CBT compared to TAU 


# Step 1: Define parameters
param <- define_parameters(
  
  age_init = 14, # assuming our cohort starts at age 14
  
  # age increases with cycles
  age = age_init + model_time,
  
  # Transition probabilities 
  p_uh_to_ps = 0.2152,                          # UHR → Psychosis
  p_uh_to_ns = 0.4465,                          # UHR → No Symptoms
  p_uh_to_ps_cbt = p_uh_to_ps*rr,
  p_ps_to_d = death_prob * 2.58 + 0.0057,       # Psychosis → Death
  p_ps_to_pp = 0.7862, #(markov paper) Psychosis → Post-Psychosis
  #p_ps_to_pp = 1-p_ps_to_d,                     
  p_pp_to_ps = 1 - exp(-(-log(1 - 0.51) / 3)),  # Recurrence over 3 years converted to annual probability
  p_pp_to_d = death_prob + 0.0057,  # Death for Post-Psychosis 
  p_pp_to_pp = 1 - (p_pp_to_ps + 0.0057 + p_pp_to_d), # Complement
  p_ns_to_d = death_prob, 
  p_ns_to_ns = 1 - p_ns_to_d,
  
  # Discount rates
  dr_costs = 0.04,       # Costs discount rate
  dr_health = 0.015,     # Health effects discount rate
  
  # Costs
  cost_uh = 4874,         # Cost per cycle in UHR
  cost_ps = 7200,         # Cost per cycle in Psychosis
  cost_pp = 5019,         # Cost per cycle in Post-Psychosis
  cost_ns = 3970,         # Cost per cycle in No Symptoms, reminder to double check this parameter
  
  #Cost TAU
  cost_tau = 1000, #NOTE TO CHANGE THIS ONE 
  
  #Cost CBT
  cost_cbt = ifelse(model_time <= 1, 1924, 1000),
  
  # Utilities
  util_uh = 0.64,        # Utility for UHR
  util_ps = 0.34,        # Utility for Psychosis
  util_pp = 0.362,       # Utility for Post-Psychosis
  util_ns = 0.8,         # Utility for No Symptoms
)

# Transition matrix for CBT
tm_cbt <- define_transition(
  C, p_uh_to_ps_cbt, p_uh_to_ns, 0, death_prob,   # UHR row
  0, 0, p_ps_to_pp, 0, C,                 # Psychosis row
  0, p_pp_to_ps, C, 0, p_pp_to_d,                 # Post-Psychosis row
  0, 0, 0, C, death_prob,                         # No Symptoms row
  0, 0, 0, 0, 1                                   # Death row
)

# Transition matrix for TAU
tm_tau <- define_transition(
  C, p_uh_to_ps, p_uh_to_ns, 0, death_prob,       # UHR row
  0, 0, p_ps_to_pp, 0, C,                 # Psychosis row
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
  init = c(1, 0, 0, 0, 0),
  cost = cost,
  effect = utility,
  method = "life-table"      # Use life-table method
)

model_tau <- run_model(
  strategy = strategy_tau,
  parameters = param,
  cycles = 10,               # Number of cycles
  init = c(1, 0, 0, 0, 0),
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
  init = c(1, 0, 0, 0, 0),
  cost = cost,
  effect = utility, 
  method = "life-table"
)

summary(results)

plot(results)

# Adding probabilistic component 

psax <- define_psa(
  cost_uh ~ gamma(mean = 200, sd = sqrt(2756)),   
  cost_ps ~ gamma(mean = 500, sd = sqrt(500)),      
  cost_pp ~ gamma(mean = 300, sd = sqrt(300)),        
  cost_ns ~ gamma(mean = 200, sd = sqrt(200)),      
  
  # Utilities 
  #util_uh ~ beta(shape1 = 142.22, shape2 = 35.56)
)


pm <- run_psa(
  model = results,
  psa = psax,
  N = 100
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

plot(pm, type = "cov")

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
