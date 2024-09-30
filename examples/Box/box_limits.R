plot(x=-1,y=-1,xlim=c(0,1),ylim=c(70,130))
add_obs(tout,obs,'black')
taus = seq(0,tout,length.out=100)[-1]
for(tau in taus){
        box = get_box_brownian_fast("BDI",theta,tout,tau,x,obs,sqrt(0.95))
        add_box(box$lower,box$upper,tout,tau,'red')
}

plot(x=-1,y=-1,xlim=c(0,1),ylim=c(70,130))
add_obs(tout,obs,'black')
taus = seq(0,tout,length.out=100)[-1]
plot_smooth_box("BDI",x,theta,obs,tout,taus,sqrt(0.95),'red')


obs = 110
plot(x=-1,y=-1,xlim=c(0,tout),ylim=c(70,130))
add_obs(tout,obs,'black')
taus = seq(0,tout,length.out=1000)[-1]
plot_smooth_box("BDI",x,theta,obs,tout,taus,sqrt(0.95),'red')
plot_smooth_box("BDI",x,theta,obs,tout,taus,sqrt(0.975),'orange')
plot_smooth_box("BDI",x,theta,obs,tout,taus,sqrt(0.99),'green')

obs = 110
plot(x=-1,y=-1,xlim=c(0,tout),ylim=c(70,130))
add_obs(tout,obs,'black')
taus = seq(0,tout,length.out=1000)[-1]
plot_smooth_box_changing_goal("BDI",x,theta,obs,tout,taus,0.9,'green')
plot_smooth_box("BDI",x,theta,obs,tout,taus,sqrt(0.9),'black')

