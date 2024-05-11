source("examples/Parameters/LV/(N=100).R")

M = 1000
print(RB("LV", x, theta, tout, lower, upper, tau, obs, M))
# 0.0003194474 tau=0.1
# 0.000322033 tau=1
# 0.00031396 tau=0.1
tau=0.1
tout

print(run_algorithm("LV", x, theta, tout, lower, upper, tau))


"LV" 
x # dataframe 
theta 
tout 
lower
upper
tau
obs
M

x
tout_list = 1:5
tout_list = 1
obs_list = get_obs("LV", x, theta, tout+1, lower, upper, tau,tout_list)
lower_list = get_lower("LV",obs_list,c(10,10))
upper_list = get_upper("LV",obs_list,c(10,10))
obs_list
lower_list
upper_list
print(RB_multi_obs("LV", x, theta, tout_list, lower_list, upper_list, tau, obs_list, M))

