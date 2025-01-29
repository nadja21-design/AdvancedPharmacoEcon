# Load packages
library(heemod)
library(dplyr)

# General death probabilities for UHR
death_prob = 0.003 # defined outside define_parameters() 
rr = 0.501 # effectiveness of CBT compared to TAU from RCT 


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
  p_ps_to_pp = 0.7862+0.2,                      # Psychosis → Post-Psychosis + remission
  p_pp_to_ps = 1 - exp(-(-log(1 - 0.51) / 3)),  # Recurrence over 3 years converted to annual probability
  p_pp_to_d = death_prob + 0.0057,              # Death for Post-Psychosis 
  p_pp_to_pp = 1 - (p_pp_to_ps + 0.0057 + p_pp_to_d), # Complement staying in post-psychosis
  p_ns_to_d = death_prob,                             # No symptoms → Death
  p_ns_to_ns = 1 - p_ns_to_d,                         # Complement for staying in no symptoms
  p_uh_to_uh = 1 - (p_uh_to_ps + p_uh_to_ns + death_prob),  # Complement for staying in UHR, MIGHT NEED TO REMOVE
   
  # Discount rates
  dr_costs = 0.04,       # Costs discount rate
  dr_health = 0.015,     # Health effects discount rate
  
  # Costs
  cost_uh = 6078,         # Cost per cycle in UHR, 3474.39
  cost_ps = 8979,         # Cost per cycle in Psychosis, 5877.82
  cost_pp = 6259,         # Cost per cycle in Post-Psychosis, 4474.39
  cost_ns = 4951,         # Cost per cycle in No Symptoms
  
  # Cost TAU
  cost_tau = 0, # reference 
  
  # Cost CBT
  cost_cbt = 2399, # should  stay consistent 
  
  # Utilities
  util_uh = 0.64,        # Utility for UHR
  util_ps = 0.34,        # Utility for Psychosis
  util_pp = 0.69,        # Utility for Post-Psychosis
  util_ns = 0.8,         # Utility for No Symptoms
)

# Step 2: define transition matrices

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

# Step 4: Define health states TAU
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
  cost_uh ~ gamma(mean = 6078, sd = 237),                 # Cost per cycle in UHR - 3.9 %
  cost_ps ~ gamma(mean = 8979, sd = 233),                 # Cost per cycle in Psychosis - 2.6%
  cost_pp ~ gamma(mean = 6259, sd = 138),                 # Cost per cycle in Post-Psychosis - 2.2%
  cost_ns ~ gamma(mean = 4951, sd = 45),                  # Cost per cycle in No Symptoms - 0.9 %
  cost_cbt ~ gamma(mean = 2399, sd = sqrt(2399)),         # Cost for CBT


  # Probabilistic utilities using beta distribution
  util_uh ~ beta(shape1 = 655, shape2 = 369),              # Utility for UHR (mean ~0.64)
  util_ps ~ beta(shape1 = 6, shape2 = 12),                 # Utility for Psychosis (mean ~0.34)
  util_pp ~ beta(shape1 = 25, shape2 = 11),                # Utility for Post-Psychosis (mean ~0.69)
  util_ns ~ beta(shape1 = 120, shape2 = 30),               # Utility for No Symptoms (mean ~0.8) 

  # Probabilistic transitions using multinomial (with explicit naming)
  p_uh_to_uh + p_uh_to_ps + p_uh_to_ns + death_prob ~ multinomial(215, 447, 335, 3),   # UHR transitions
  p_ps_to_d ~ binomial(prob = 0.013, size = 1000),                                     # Psychosis transitions
  p_pp_to_pp + p_pp_to_ps + p_pp_to_d ~ multinomial(510, 362, 5),                      # Post-Psychosis transitions
  p_ns_to_d ~ binomial(prob = 0.001, size = 1000),                                      # No Symptoms transitions
  
  rr ~ lognormal(meanlog = log(0.501), sdlog = 0.017)
  )

set.seed(123) 

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

plot(pm, type = "ce")

# Cost-effectiveness acceptability curves or EVPI:

plot(pm, type = "ac", max_wtp = 10000, log_scale = FALSE)
plot(pm, type = "evpi", max_wtp = 10000, log_scale = FALSE)

# A covariance analysis can be performed on strategy results:

plot(pm, type = "cov") 
plot(pm, type = "ac", log_scale = FALSE, max_wtp = 10000)

# Or on the difference between strategies:

plot(pm, type = "cov", diff = TRUE, threshold = 5000)

# ggplot2 syntax with WTP threshold 

library(ggplot2)

wtp_line = 20000
  
plot(pm, type = "ce") +
  xlab("QALYs gained") +
  ylab("Additional cost") +
  geom_abline(slope = wtp_line, intercept = 0, color = "black", linetype="dashed")+
  scale_color_brewer(
    name = "Strategy",
    palette = "Set1"
  ) +
  theme_minimal()

#####