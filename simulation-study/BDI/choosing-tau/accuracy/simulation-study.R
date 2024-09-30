### PVR = 0.90
pdf("simulation-study/BDI/choosing-tau/accuracy/equilibrium.pdf")
source('examples/Parameters/BDI/Equilibrium (N=100).R')
tau=0.05
plot_PVR_target_goal_tau("BDI",x,theta,obs_list[1,1],tout_list[1])
plot_PVR_target_goal_width("BDI",x,theta,obs_list[1,1],tout_list[1],tau)
plot_PVR_target_goal("BDI",x,theta,obs_list[1,1],tout_list[1])
dev.off()

### outlier, goal = 0.90
pdf("simulation-study/BDI/choosing-tau/accuracy/equilibrium-outlier.pdf")
source("examples/Parameters/BDI/Equilibrium (outlier 2SD).R")
tau=0.1
plot_PVR_target_goal_tau("BDI",x,theta,obs_list[1,1],tout_list[1])
plot_PVR_target_goal_width("BDI",x,theta,obs_list[1,1],tout_list[1],tau)
plot_PVR_target_goal("BDI",x,theta,obs_list[1,1],tout_list[1])
dev.off()

source('examples/Parameters/BDI/Equilibrium (N=100).R')
source('examples/Parameters/BDI/Drift/Drift (N=100).R')

iterations = 30000
num_particles = 65

tout_list = c(0.25,0.5,0.75,1)

tout_list = seq(1,10)
tout_list = seq(0,30,length.out=11)[-1]
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[,1] = c(96, 102, 113, 106, 93, 102, 81, 83, 90, 89)
# 1.5, 53 44 42
z1.start=Sys.time()
z1 = mcmc_RB("BDI",x,theta,obs_list,tout_list,20000,30)
z1.end=Sys.time()
z2.start=Sys.time()
z2 = mcmc_frac("BDI",x,theta,obs_list,tout_list,iterations,700)
z2.end=Sys.time()
c(135,130,100)
c(154,141,116)
c(294,318,281)
z1.end-z1.start
z2.end-z2.start

t1 = difftime(z1.end,z1.start,units='mins')[[1]]
t2 = difftime(z2.end,z2.start,units='mins')[[1]]

z1$acceptance_rate
z2$acceptance_rate
var(z1$ll)
var(z2$ll)
plot(z1$log_posterior[,1],type='l')
plot(z2$log_posterior[,1],type='l')

signif((effectiveSize(z1$log_posterior)/t1),3)
signif((effectiveSize(z2$log_posterior)/t2),3)
signif((effectiveSize(z1$log_posterior)/t1)/(effectiveSize(z2$log_posterior)/t2),3)

colnames(z1$log_posterior) = c("Birth_rate","Death_rate","Immigration")
colnames(z2$log_posterior) = c("Birth_rate","Death_rate","Immigration")
pdf("simulation-study/BDI/multi-obs/obs=30/plot.pdf")
pairs(exp(z1$log_posterior),cex=0.5,cex.axis=1.5,cex.lab=1.5,main="RB",col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),pch=19)
pairs(exp(z2$log_posterior),cex=0.5,cex.axis=1.5,cex.lab=1.5,main="Trivial",col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),pch=19)
dev.off()
library(GGally)
ggpairs(exp(z1$log_posterior))+ theme(text = element_text(size = 30))


apply(exp(z1$log_posterior),mean,MARGIN = 2)
apply(exp(z2$log_posterior),mean,MARGIN = 2)
get_mu("BDI",x,theta)
get_mu("BDI",x,apply(exp(z1$log_posterior),mean,MARGIN = 2))
get_mu("BDI",x,apply(exp(z2$log_posterior),mean,MARGIN = 2))


library(rgl)
plot3d(exp(z1$log_posterior),size=10)

mean(exp(z1$ll))
mean(exp(z2$ll))

hist(exp(z1$log_posterior[,1]),breaks=20)

plot(exp(z1$log_posterior[,1]),exp(z1$log_posterior[,2]))
plot(exp(z1$log_posterior[,1]),exp(z1$log_posterior[,3]))
plot(exp(z1$log_posterior[,2]),exp(z1$log_posterior[,3]))

mean(exp(z1$log_posterior[,1]))
mean(exp(z1$log_posterior[,2]))
mean(exp(z1$log_posterior[,3]))

var(z1$ll)
cov(z1$log_posterior)
z1$acceptance_rate

head(as.data.frame(exp(z1$log_posterior)))

#####
ggplot(as.data.frame(exp(z1$log_posterior)))+geom_density2d_filled(aes(x=Birth_rate,y=Death_rate))
ggplot(as.data.frame(exp(z2$log_posterior)))+geom_density2d_filled(aes(x=Birth_rate,y=Death_rate))

par(list(mfrow=c(1,3)))
par(mar=c(2,2,1,1))

pdf("simulation-study/BDI/multi-obs/obs=30/posterior.pdf")
plot(-1,xlim=c(-6,6),ylim=c(0,1.2),xlab='(log-scale)',ylab='density',main='Posterior density estimate',cex.main=2,cex.axis=2,cex.lab=2)
lines(density(z1$log_posterior[,1]),col='blue',lwd=2)
lines(density(z2$log_posterior[,1]),col='blue',lty='dashed',lwd=2)
lines(density(z1$log_posterior[,2]),col='red',lwd=2)
lines(density(z2$log_posterior[,2]),col='red',lty='dashed',lwd=2)
lines(density(z1$log_posterior[,3]),col='green',lwd=2)
lines(density(z2$log_posterior[,3]),col='green',lty='dashed',lwd=2)
legend(x=0.5,y=1.1,c("Birth rate (RB)","Death rate (RB)","Immigration rate (RB)"),col=c('blue','red','green'),lty=c('solid','solid','solid'),lwd=2)
legend(x=0.5,y=0.8,c("Birth rate (Trival)","Death rate (Trival)","Immigration rate (Trival)"),col=c('blue','red','green'),lty=c('dashed','dashed','dashed'),lwd=2)
dev.off()

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("blue","black", "green", "orange", "pink")
legend(x = "top",inset = 0,
        legend = c("Fabricated Metal", "Iron and Steel", "Paper","Beverages", "Tobacco"), 
        col=plot_colors, lwd=5, cex=.5, horiz = TRUE)


Birth_rate = c(effectiveSize(z1$log_posterior)[1],effectiveSize(z2$log_posterior)[1])
Death_rate = c(effectiveSize(z1$log_posterior)[2],effectiveSize(z2$log_posterior)[2])
Immigration_rate = c(effectiveSize(z1$log_posterior)[3],effectiveSize(z2$log_posterior)[3])

df = data.frame(Iterations=c(20000,20000),'Birth rate'=Birth_rate,'Death rate'=Death_rate,'Immigration rate'=Immigration_rate)
xtable(df)

addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- paste0(paste0('&& \\multicolumn{2}{c}{', "ESS", '}', collapse=''), '\\\\')

print(xtable(df),add.to.row=addtorow)


particles = seq(20,40,by=2)
results = data.frame(matrix(nrow=length(particles),ncol=6))

names(results) = c("Particles","Var_Logpost","Accept","ESS_birth","ESS_death","ESS_immigration")
for(i in 1:length(particles)){
        z1.start=Sys.time()
        z1 = mcmc_RB("BDI",x,theta,obs_list,tout_list,20000,particles[i])
        z1.end=Sys.time()
        t1 = difftime(z1.end,z1.start,units='mins')[[1]]
        ESS = signif((effectiveSize(z1$log_posterior)/t1),3)
        results[i,] = c(particles[i],var(z1$ll),z1$acceptance_rate,ESS)
}
results
plot(results$Particles,results$ESS_birth,type='l')

results_new = rbind(results,results2)
results_new = results_new[order(results_new$Particles),]
plot(z1$log_posterior[,1],type='l')
goals = 1-1/(10*2^(1:5))
# Create sequence of goals
for(i in 1:length(goals)){
        tau_list = get_tau_list("BDI",x,theta,obs_list,tout_list,goals[i])
        box_list = get_box_brownian_fast_list("BDI",theta,tout_list,tau_list,x,obs_list,goals[i])
}
plot((box_list$upper-box_list$lower),)
# Get box
# Get width
# Plot (widths, goals)

tau_list = get_tau_list("BDI",x,theta,obs_list,tout_list,sqrt(0.95))
box = get_box_brownian_fast_list("BDI",theta,tout_list,tau_list,x,obs_list,sqrt(0.95))

B = 5000
a_results = rep(0,B)
b_results = rep(0,B)
c_results = rep(0,B)

a_seq = seq(0,0.6,length.out=B)
b_seq = seq(0,0.6,length.out=B)
c_seq = seq(0,20,length.out=B)

for(i in 1:B){
        a_results[i] = RB_list("BDI",x,c(a_seq[i],theta[2],theta[3]),tout_list,box$lower,box$upper,tau_list,obs_list,30) 
}
for(i in 1:B){
        b_results[i] = RB_list("BDI",x,c(theta[1],b_seq[i],theta[3]),tout_list,box$lower,box$upper,tau_list,obs_list,30) 
}
for(i in 1:B){
        c_results[i] = RB_list("BDI",x,c(theta[1],theta[2],c_seq[i]),tout_list,box$lower,box$upper,tau_list,obs_list,30) 
}
pdf('tuning.pdf')
plot(a_seq,a_results,xlab='Birth rate',ylab='Marginal likelihood (30 particles)',cex.lab=2)
plot(b_seq,b_results,xlab='Death rate',ylab='Marginal likelihood (30 particles)',cex.lab=2)
plot(c_seq,c_results,xlab='Immigration rate',ylab='Marginal likelihood (30 particles)',cex.lab=2)
dev.off()
var(z1$ll)

d_results = rep(0,20000)
for(i in 1:length(d_results)){
        print(i)
        d_results[i] = RB_list("BDI",x,theta,tout_list,box$lower,box$upper,tau_list,obs_list,100) 
}
var(d_results[d_results>-1e+20])
var(exp(d_results))/mean(exp(d_results))^2
hist(d_results)
exp(1)-1
hist(exp(d_results))
head(rev(sort(d_results)),20)
#############################

prev_results = results

particles = seq(60,10,by=-10)
results = data.frame(matrix(nrow=length(particles),ncol=5))
results = data.frame(matrix(nrow=length(particles),ncol=2))
names(results) = c("Particles","Var_Logpost","Accept","ESS_birth","ESS_death","ESS_immigration")
tau_list = get_tau_list("BDI",x,theta,obs_list,tout_list,sqrt(0.95))
box_list = get_box_brownian_fast_list("BDI",theta,tout_list,tau_list,x,obs_list,sqrt(0.95))
for(i in 1:length(particles)){
        d_results = rep(0,5000)
        print(i)
        for(j in 1:length(d_results)){
                d_results[j] = RB_list("BDI",x,theta,tout_list,box_list$lower,box_list$upper,tau_list,obs_list,particles[i]) 
        }
        results[i,] = c(particles[i],var(d_results[d_results>-1e+20]))
}

for(i in 1:length(particles)){
        z1.start=Sys.time()
        z1 = mcmc_RB("BDI",x,theta,obs_list,tout_list,200000,particles[i])
        z1.end=Sys.time()
        t1 = difftime(z1.end,z1.start,units='mins')[[1]]
        ESS = signif((effectiveSize(z1$log_posterior)/t1),3)
        results[i,] = c(particles[i],z1$acceptance_rate,ESS)
}
plot(results,type='o')

cbind(prev_results$X2,results)

colnames(results) = c('particles','acceptance',"ESS1",'ESS2','ESS3')
colnames(prev_results) = c('particles','var')


###############
E_cpu = log(mean(exp(get_q(100,1))))
r = mean(is.infinite(get_q(1000,1))) # prob of particle missing at least 1 obs
s = mean(is.infinite(get_q(1000,10)))
q = exp(log(1-exp(log(1-r)/(nout-1)))/1)
h = exp(log(1-exp(log(1-s)/(nout-1)))/20)

1-(1-q)^(nout-1)
1-(1-h^10)^(nout-1)

q
h


