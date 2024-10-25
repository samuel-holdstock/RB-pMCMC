source("examples/Parameters/SEIR_N/Equilibrium (N=1000).R")
lower = c(1,0,0)
upper = c(5,5,80)

obs_index = state_to_index(obs,lower,upper)
tau = 20

data = sim_data("SEIR_N",x,theta, tout, lower, upper, tau)$data
data_chain = cbind(data,N - apply(data[,-1],1,sum))
plot_chain(data_chain,c("Susceptible","Exposed","Infected","Removed"),c("blue","orange","darkgreen","red"),ylim=c(0,1000))
add_box(lower,upper,tout,tau,c("blue","orange","darkgreen"))
add_obs(tout,obs,c("blue","orange","darkgreen"))


state_to_index(lower+c(0,0,1),lower,upper)
index_to_state(10209,lower,upper)



noise_lower = c(5,10,15)
x = c(10,15,20)
noise_upper = c(15,30,25)

lower = c(1,2,3)
upper = c(20,40,30)

get_noise_index = function(noise_lower,noise_upper,lower,upper,coords){
  noise_lower_index = state_to_index(noise_lower,lower,upper)
  num_params = length(x)
  index_jumps = cumprod(rev(upper-lower)+1)
  index_jumps = c(1,index_jumps[-num_params])
  index = noise_lower_index + sum(index_jumps*coords)
  print(index_to_state(index,lower,upper))
}
get_noise_index(noise_lower,noise_upper,lower,upper,c(0,0,0))
index_to_state(get_noise_indices(noise_lower,noise_upper,lower,upper,c(2,3,1)),lower,upper)
get_noise_indices(get_noise_indices(noise_lower,noise_upper,lower,upper,c(2,3,1)),lower,upper)

get_noise_indices(noise_lower,noise_upper,lower,upper)

noise_upper=c(10,15,20)
noise_upper-noise_lower
