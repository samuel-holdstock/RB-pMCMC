get_obs = function(model_name,start,theta,times){
    num_obs = length(times)
    data = sim_data_frac(model_name,x,theta,times[num_obs])$data
    obs_list = data.frame(matrix(ncol=length(start)))
    num_data = nrow(data)
    #times = c(0,times)
    count = 1
    for(i in 1:num_data){
        for(j in count:num_obs){
            next_time = times[j]
            if(next_time>=data[i,1] && (i==num_data || data[i+1,1]>next_time)){
                x_row = data[i,]
                x_row[1] = next_time
                obs_list[j,] = x_row[-1]
                count = count + 1
            }
        }
    }    
    return(as.matrix(unname(obs_list)))
}

get_lower = function(model_name,obs_list,widths){
    num_obs = nrow(obs_list)
    dim = ncol(obs_list)
    lower_list = data.frame(matrix(ncol=dim))
    for(i in 1:num_obs){
        lower_vec = rep(0,dim)
        for(j in 1:dim){
            lower_vec[j] = floor(obs_list[i,][j]-widths[j]/2)
        }        
        lower_list[i,] = lower_vec
    }
    return(as.matrix(unname(lower_list)))
}

get_upper = function(model_name,obs_list,widths){
    num_obs = nrow(obs_list)
    dim = ncol(obs_list)
    upper_list = data.frame(matrix(ncol=dim))
    for(i in 1:num_obs){
        upper_vec = rep(0,dim)
        for(j in 1:dim){
            upper_vec[j] = ceiling(obs_list[i,][j]+widths[j]/2)
        }        
        upper_list[i,] = upper_vec
    }
    return(as.matrix(unname(upper_list)))
}

get_box = function(model_name,theta,tout_list,tau,x0,obs_list,pvr_goal){
    num_obs = nrow(obs_list)
    dim = ncol(obs_list)
    lower_list = data.frame(matrix(ncol=dim))
    upper_list = data.frame(matrix(ncol=dim))
    x = x0
    prev_time = 0
    for(i in 1:num_obs){
        box = get_box_brownian_fast(model_name,theta,tout_list[i]-prev_time,tau,x,obs_list[i,],pvr_goal)
        lower_list[i,] = box$lower
        upper_list[i,] = box$upper
        x = obs_list[i,]
        prev_time = tout_list[i]
    }
    return(list(lower=as.matrix(unname(lower_list)), upper=as.matrix(unname(upper_list))))
}

plot_RB_MJP = function(model_name, x0, theta, tout_list,lower_list,upper_list, tau, obs_list, M, species_names,species_col,ylims){
    num_obs = nrow(obs_list)
    num_species = ncol(obs_list)
    tout = tout_list[num_obs]
    plot(1,xlim=c(0,tout),ylim=ylims,main="Reaction network",xlab="Time",ylab="Species count")
    for(i in 1:num_obs){
        lower = lower_list[i,]
        upper = upper_list[i,]
        obs = obs_list[i,]
        tout = tout_list[i]
        add_box(lower,upper,tout,tau,species_col)
        add_obs(tout,obs,species_col)
    }

    for(i in 1:M){
        x = x0
        prev_time = 0
        for(j in 1:(num_obs)){
            if(j>1){
                x = obs_list[j-1,]
            }
            tout = tout_list[j]-prev_time
            data = sim_data_frac(model_name,x,theta,tout)$data
            data_chain = get_chain(data,tout)
            data_chain[,1] = data_chain[,1]+prev_time
            add_chain(data_chain,species_names,alpha(species_col,0.3))
            prev_time = tout_list[j]
        }
    }
    for(i in 1:num_species){
        prev_time = 0
        x = x0
        for(j in 1:num_obs){
            x_coords=c(prev_time,tout_list[j])
            mu = get_mu(model_name,x,theta)
            forecast = x[i] + mu[i]*(tout_list[j]-prev_time)
            y_coords=c(x[i],forecast)
            prev_time = tout_list[j]
            x = obs_list[j,]
            lines(x_coords,y_coords,lty='dashed',col=species_col[i])
        }
    }
    legend("topleft",legend=species_names,col=species_col,inset=0.05,lty=rep(1,length(species_names)))
}

get_study_rb = function(model_name,x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau){
    rb_sims = list(sims=rep(0,N),times=rep(0,N))
    for(i in 1:N){
        start.time = Sys.time()
        rb_sims$sims[i] = RB_list(model_name,x, theta, tout_list, lower_list, upper_list, tau,obs_list,M)
        end.time = Sys.time()
        rb_sims$times[i] = end.time-start.time
    }
    return(rb_sims)
}

get_study_frac = function(model_name,x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau){
    frac_sims = list(sims=rep(0,N),times=rep(0,N))
    for(i in 1:N){
        start.time = Sys.time()
        frac_sims$sims[i] = frac_list(model_name,x, theta, tout_list, obs_list,M)
        end.time = Sys.time()
        frac_sims$times[i] = end.time-start.time
    }
    return(frac_sims)
}

plot_study_estimates_tau = function(rb,frac,barwidth,x_max=1,tau){
    rb_est = rb$sims
    frac_est = frac$sims

    breaks = (seq(0,max(c(rb_est,frac_est))+2*barwidth,by=barwidth)-barwidth/2)
    p1 = hist(rb_est,breaks=breaks,plot=F)
    p1$counts=p1$counts/sum(p1$counts)
    p2 = hist(frac_est,breaks=breaks,plot=F)
    p2$counts=p2$counts/sum(p2$counts)

    q1 = p1$breaks[which(cumsum(p1$counts)>=x_max)[1]+1]
    q2 = p2$breaks[which(cumsum(p2$counts)>=x_max)[1]+1]
    xlim = c(-barwidth/2,max(c(q1,q2)))

    plot(p1,main='Estimates of likelihood',xlab='Likelihood',ylab='Probability Distribution',xlim=xlim,ylim=c(0,1),
    cex.lab=1.5,cex.main=1.5,cex.axis=1.5,yaxt='n',col=rgb(1,0,0,1/4))
    plot(p2,xlim=xlim,ylim=c(0,1),add=T,col=rgb(0,0,1,1/4))
    axis(2,cex.axis=1.5)

    abline(v=mean(rb_est),col=rgb(1,0,0,1/4),lty='dashed')
    abline(v=mean(frac_est),col=rgb(0,0,1,1/4),lty='dashed')

    times_rb = rb$times
    times_frac = frac$times
    rb_eff = formatC(var(rb_est)*sum(times_rb),format='e')
    frac_eff = formatC(var(frac_est)*sum(times_frac),format='e')
    eff_ratio = round((var(frac_est)*sum(times_frac))/(var(rb_est)*sum(times_rb)),2)

    mtext(side = 3, line = 0.25, adj = 0.5, paste("tau =",tau),cex = 1.5)
    mtext(side = 3, line = -10, adj = 0.8, paste("RB Inefficiency =",rb_eff),cex = 1,col='red')
    mtext(side = 3, line = -11, adj = 0.8, paste("Frac Inefficiency =",frac_eff),cex = 1,col='blue')
    mtext(side = 3, line = -13, adj = 0.8, paste("Ratio =",eff_ratio),cex = 1.5,col='black')
    legend(x=xlim[2]*0.7,y=0.25,c("RB","Fraction"),c("red","blue"))    
}

plot_study_estimates_M = function(rb,frac,barwidth,x_max=1){
    rb_est = rb$sims
    frac_est = frac$sims

    breaks = (seq(0,max(c(rb_est,frac_est))+2*barwidth,by=barwidth)-barwidth/2)

    p1 = hist(rb_est,breaks=breaks,plot=F)
    p1$counts=p1$counts/sum(p1$counts)
    p2 = hist(frac_est,breaks=breaks,plot=F)
    p2$counts=p2$counts/sum(p2$counts)

    q1 = p1$breaks[which(cumsum(p1$counts)>=x_max)[1]+1]
    q2 = p2$breaks[which(cumsum(p2$counts)>=x_max)[1]+1]
    xlim = c(-barwidth/2,max(c(q1,q2)))

    plot(p1,main='Estimates of likelihood',xlab='Likelihood',ylab='Probability Distribution',xlim=xlim,ylim=c(0,1),
    cex.lab=1.5,cex.main=1.5,cex.axis=1.5,yaxt='n',col=rgb(1,0,0,1/4))
    plot(p2,xlim=xlim,ylim=c(0,1),add=T,col=rgb(0,0,1,1/4))
    axis(2,cex.axis=1.5)

    abline(v=mean(rb_est),col=rgb(1,0,0,1/4),lty='dashed')
    abline(v=mean(frac_est),col=rgb(0,0,1,1/4),lty='dashed')

    times_rb = rb$times
    times_frac = frac$times
    rb_eff = formatC(var(rb_est)*sum(times_rb),format='e')
    frac_eff = formatC(var(frac_est)*sum(times_frac),format='e')
    eff_ratio = round((var(frac_est)*sum(times_frac))/(var(rb_est)*sum(times_rb)),2)

    mtext(side = 3, line = 0.25, adj = 0.5, paste("M =",M),cex = 1.5)
    mtext(side = 3, line = -10, adj = 0.8, paste("RB Inefficiency =",rb_eff),cex = 1,col='red')
    mtext(side = 3, line = -11, adj = 0.8, paste("Frac Inefficiency =",frac_eff),cex = 1,col='blue')
    mtext(side = 3, line = -13, adj = 0.8, paste("Ratio =",eff_ratio),cex = 1.5,col='black')
    legend(x=xlim[2]*0.7,y=0.25,c("RB","Fraction"),c("red","blue"))
}

plot_study_estimates_goal = function(rb,frac,barwidth,x_max=1,goal){
    rb_est = rb$sims
    frac_est = frac$sims

    breaks = (seq(0,max(c(rb_est,frac_est))+2*barwidth,by=barwidth)-barwidth/2)
    p1 = hist(rb_est,breaks=breaks,plot=F)
    p1$counts=p1$counts/sum(p1$counts)
    p2 = hist(frac_est,breaks=breaks,plot=F)
    p2$counts=p2$counts/sum(p2$counts)

    q1 = p1$breaks[which(cumsum(p1$counts)>=x_max)[1]+1]
    q2 = p2$breaks[which(cumsum(p2$counts)>=x_max)[1]+1]
    xlim = c(-barwidth/2,max(c(q1,q2)))

    plot(p1,main='Estimates of likelihood',xlab='Likelihood',ylab='Probability Distribution',xlim=xlim,ylim=c(0,1),
    cex.lab=1.5,cex.main=1.5,cex.axis=1.5,yaxt='n',col=rgb(1,0,0,1/4))
    plot(p2,xlim=xlim,ylim=c(0,1),add=T,col=rgb(0,0,1,1/4))
    axis(2,cex.axis=1.5)

    abline(v=mean(rb_est),col=rgb(1,0,0,1/4),lty='dashed')
    abline(v=mean(frac_est),col=rgb(0,0,1,1/4),lty='dashed')

    times_rb = rb$times
    times_frac = frac$times
    rb_eff = formatC(var(rb_est)*sum(times_rb),format='e')
    frac_eff = formatC(var(frac_est)*sum(times_frac),format='e')
    eff_ratio = round((var(frac_est)*sum(times_frac))/(var(rb_est)*sum(times_rb)),2)
    achieved = round(1-var(rb_est)/var(frac_est),2)

    mtext(side = 3, line = 0.25, adj = 0.5, paste("goal =",goal),cex = 1.5)
    mtext(side = 3, line = -1, adj = 0.5, paste("achieved =",achieved),cex = 1.5)
    mtext(side = 3, line = -4, adj = 0.8, paste("RB Inefficiency =",rb_eff),cex = 1,col='red')
    mtext(side = 3, line = -5, adj = 0.8, paste("Frac Inefficiency =",frac_eff),cex = 1,col='blue')
    mtext(side = 3, line = -7, adj = 0.8, paste("Ratio =",eff_ratio),cex = 1.5,col='black')
    legend(x=xlim[2]*0.7,y=0.55,c("RB","Fraction"),c("red","blue"))    
}

plot_study_estimates_generic = function(rb,frac,barwidth,x_max=1,subtitle){
    rb$sims = exp(rb$sims)
    frac$sims = exp(frac$sims)

    rb_est = rb$sims
    frac_est = frac$sims

    breaks = (seq(0,(max(c(rb_est,frac_est))+2*barwidth),by=barwidth)-barwidth/2)
    p1 = hist(rb_est,breaks=breaks,plot=F)
    p1$counts=p1$counts/sum(p1$counts)
    p2 = hist(frac_est,breaks=breaks,plot=F)
    p2$counts=p2$counts/sum(p2$counts)

    q1 = p1$breaks[which(cumsum(p1$counts)>=x_max)[1]+1]
    q2 = p2$breaks[which(cumsum(p2$counts)>=x_max)[1]+1]
    xlim = c(-barwidth/2,max(c(q1,q2)))

    plot(p1,main='Estimates of likelihood',xlab='Likelihood',ylab='Probability Distribution',xlim=xlim,ylim=c(0,1),
    cex.lab=1.5,cex.main=1.5,cex.axis=1.5,yaxt='n',col=rgb(1,0,0,1/4))
    plot(p2,xlim=xlim,ylim=c(0,1),add=T,col=rgb(0,0,1,1/4))
    axis(2,cex.axis=1.5)

    abline(v=mean(rb_est),col=rgb(1,0,0,1/4),lty='dashed')
    abline(v=mean(frac_est),col=rgb(0,0,1,1/4),lty='dashed')

    times_rb = rb$times
    times_frac = frac$times
    rb_eff = formatC(var(rb_est)*sum(times_rb),format='e')
    frac_eff = formatC(var(frac_est)*sum(times_frac),format='e')
    eff_ratio = round((var(frac_est)*sum(times_frac))/(var(rb_est)*sum(times_rb)),2)

    mtext(side = 3, line = 0.25, adj = 0.5, subtitle,cex = 1.5)
    # mtext(side = 3, line = -4, adj = 0.8, paste("RB Inefficiency =",rb_eff),cex = 1,col='red')
    # mtext(side = 3, line = -5, adj = 0.8, paste("Frac Inefficiency =",frac_eff),cex = 1,col='blue')
    # mtext(side = 3, line = -7, adj = 0.8, paste("Ratio =",eff_ratio),cex = 1.5,col='black')
    legend(x=xlim[2]*0.7,y=0.55,c("RB","Fraction"),c("red","blue"))    
}

get_row_names = function(){
    row_names = c('Probability', 'RB variance', 'Trivial variance','PVR tout', 'PVR tau')
    return(row_names)
}

plot_study_estimates_M_time = function(rb_sims,frac_sims,num_quants=10){
    rb_times = rb_sims$times
    frac_times = frac_sims$times
    
    quants = seq(0,0.95,length.out=num_quants)
    ymin = min(quantile(rb_times,min(quants)),quantile(frac_times,min(quants)))
    ymax = max(quantile(rb_times,max(quants)),quantile(frac_times,max(quants)))

    plot(quantile(rb_times,quants),quants, col='red',type='l',xlim=c(ymin,ymax),main='ECDF of time per iteration',ylab='Quantiles (0-0.95)',xlab='Time',
    cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
    lines(quantile(frac_times,quants),quants,type='l',col='blue')
    abline(v=mean(rb_times),lty='dashed',col='red')
    abline(v=mean(frac_times),lty='dashed',col='blue')
}

exp_Q_obs = function(model_name,lower,upper,theta,tau,obs){
    Q = get_coffin_matrix(model_name,lower,upper,theta)
    v = matrix(0,ncol=ncol(Q))
    obs_index = state_to_index(obs,lower,upper)+1
    v[obs_index] = 1
    vTexpQ = vT_exp_Q(v,t(tau*Q),1e-40,F,T,F)
    return(vTexpQ)
}

exp_Q_x = function(model_name,lower,upper,theta,tau,x){
    Q = get_coffin_matrix(model_name,lower,upper,theta)
    v = matrix(0,ncol=ncol(Q))
    x_index = state_to_index(x,lower,upper)+1
    v[x_index] = 1
    vexpQ = v_exp_Q(v,(tau*Q),1e-40,F,T,F)
    return(vexpQ)
}

get_PVR = function(model_name,lower_list,upper_list,theta,tau,tout_list,x0,obs_list,lower_limit_list,upper_limit_list,M){
    num_obs = length(tout_list)
    x = x0
    col_names = c('probability','rb_variance','frac_variance','PVR_tout','PVR_tau')
    df = data.frame(matrix(ncol=length(col_names),nrow=0))
    colnames(df) = col_names
    for(i in num_obs){
        lower = lower_list[i,]
        upper = upper_list[i,]
        tout = tout_list[i]
        obs = obs_list[i,]
        lower_limit = lower_limit_list[i,]
        upper_limit = upper_limit_list[i,]

        ai = as.numeric(exp_Q_x(model_name,lower_limit,upper_limit,theta,tout-tau,x)) # a_i
        bi = as.numeric(exp_Q_obs(model_name,lower_limit,upper_limit,theta,tau,obs)) # b_i
        qi = as.numeric(exp_Q_obs(model_name,lower,upper,theta,tau,obs)) # b_i - delta_i
        
        ai = ai[-length(ai)]
        bi = bi[-length(bi)]
        qi = qi[-length(qi)]
        
        p = sum(ai*bi)
        max_reduction_tau = sum(ai*bi*(1-bi))

        lower_index = state_to_index(lower,lower_limit,upper_limit)+1
        upper_index = state_to_index(upper,lower_limit,upper_limit)+1
        ai = ai[lower_index:upper_index]
        bi = bi[lower_index:upper_index]
        deltai = bi-qi

        probability = p
        frac_variance = p*(1-p)
        rb_variance = p*(1-p) - sum(ai*bi*(1-bi)) + sum(ai*deltai*(1-deltai))
        rb_reduction = sum(ai*bi*(1-bi)) - sum(ai*deltai*(1-deltai))
        PVR_tout = rb_reduction/frac_variance
        PVR_tau = rb_reduction/max_reduction_tau
        df[i,] = c(probability,rb_variance/M,frac_variance/M,PVR_tout,PVR_tau)
        x = obs_list[i,]
    }
    return(df)
}

plot_table = function(title,row_names,row_data,...){
    plot.new()
    legend(
    ncol = 2L, title = title,
    legend = c(
        row_names,
        row_data
    ),...
    )
}

add_table = function(title,row_names,row_data,...){
    legend(
    ncol = 2L, title = title,
    legend = c(
        row_names,
        row_data
    ),...
    )
}

get_row_data_a = function(row_data_a){
     row_data_a = signif(row_data_a,5)
     row_data_a$rb_variance = formatC(row_data_a$rb_variance,format='e')
     row_data_a$frac_variance = formatC(row_data_a$frac_variance,format='e')
     row_data_a = c(row_data_a)
     return(row_data_a)
}

get_row_data_b = function(rb_sims,frac_sims,M,N,row_data_a){
    row_data_b = c(
        M,
        N,
        signif(mean(rb_sims$sims),5),
        signif(mean(frac_sims$sims),5),
        formatC(var(rb_sims$sims),format='e'),
        formatC(var(frac_sims$sims),format='e'),
        signif(mean(rb_sims$times),5),
        signif(mean(frac_sims$times),5),
        formatC(1/mean(1/var(rb_sims$sims))*mean(rb_sims$times),format='e'),
        formatC(1/mean(1/var(frac_sims$sims))*mean(frac_sims$times),format='e'),
        signif(1/mean(1/var(frac_sims$sims))*mean(frac_sims$times)/(1/mean(1/var(rb_sims$sims))*mean(rb_sims$times)),5))
    return(row_data_b)
}

get_row_names = function(){
    row_names_a = c('Probability', 'RB variance', 'Naive variance','PVR tout', 'PVR tau')
    row_names_b = c('Particles','Iterations','RB probability','Trivial probability','RB variance','Trivial variance','RB time per iteration', 'Naive time per iteration','RB inefficiency','Naive inefficiency','Efficiency ratio')
    return(list(a=row_names_a,b=row_names_b))
}

get_row_data = function(model_name,lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims){
    row_data_a = get_PVR(model_name,lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M)[1,]
    row_data_b = get_row_data_b(rb_sims,frac_sims,M,N,row_data_a)
    row_data_a = get_row_data_a(row_data_a)
    return(list(a=row_data_a,b=row_data_b))
}

plot_tables = function(model_name,lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims){
    rb_sims$sims = exp(rb_sims$sims)
    frac_sims$sims = exp(frac_sims$sims)

    row_names = get_row_names()
    row_names_a = row_names$a
    row_names_b = row_names$b

    row_data = get_row_data(model_name,lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
    row_data_a = row_data$a
    row_data_b = row_data$b

    plot_table("Results via matrix exponential",row_names_a,row_data_a,x='top',cex=1.2,text.width=strwidth(toString(rep('1',14))),x.intersp=6)
    add_table("Results via simulation",row_names_b,row_data_b,x=-0.3,y=0.65,cex=1.2,text.width=strwidth(toString(rep('1',12))),x.intersp=8.1)    
}

# get_var_big_box_tau = function(model_name,x0,theta,xt,tout,tau){
#   mu = get_mu(model_name,x0,theta)
#   sig2 = get_var(model_name,x0,theta)
#   RB_variance = 1/sqrt(4*pi*sig2*tau)*1/sqrt(2*pi*sig2*(tout-tau/2))*exp(-1/2*(xt-(x0+mu*tout))^2/(sig2*(tout-tau/2)))-
#     1/(2*pi*sig2*tout)*exp(-1/2*(xt-x0-mu*tout)^2/(sig2*tout/2))
#   return(RB_variance)
#   # p = get_probability_hit(x,obs,tout,0,sqrt(60))
#   # Q1 = get_Q1_brownian_fast(0, theta, tout, 0.1, x, obs, 0, 60)
#   # return(p*(1-p)-Q1)
# }

# get_PVR_big_box_tau = function(model_name,x0,theta,xt,tout,tau){
#   mu = get_mu(model_name,x0,theta)
#   sig2 = get_var(model_name,x0,theta)
#   RB_variance = 1/sqrt(4*pi*sig2*tau)*1/sqrt(2*pi*sig2*(tout-tau/2))*exp(-1/2*(xt-(x0+mu*tout))^2/(sig2*(tout-tau/2)))-
#   1/(2*pi*sig2*tout)*exp(-1/2*(xt-x0-mu*tout)^2/(sig2*tout/2))
#   p = 1/sqrt(2*pi*sig2*tout)*exp(-1/2*(xt-x0-mu*tout)^2/(sig2*tout))
#   PVR = (p*(1-p)-RB_variance)/(p*(1-p))
#   return(PVR)
#   # p = get_probability_hit(x,obs,tout,0,sqrt(60))
#   # Q1 = get_Q1_brownian_fast(0, theta, tout, 0.1, x, obs, 0, 60)
#   # return(Q1/(p*(1-p)))
# }

plot_var_big_box_tau = function(model_name,x,theta,obs,tout){
    taus = seq(0,tout,length.out=1000)[-1]
    best_variance = rep(0,length(taus))
    for(i in 1:length(taus)){
        tau = taus[i]
        best_variance[i] = get_var_big_box_tau(model_name,x,theta,obs,tout,tau)
    }
    plot(taus,best_variance,type='l',xlab='Tau',ylab='Infinite box var',cex=1.5,cex.axis=1.5,cex.lab=1.5)
}

plot_PVR_big_box_tau = function(model_name,x,theta,obs,tout){
    taus = seq(0,tout,length.out=1000)[-1]
    best_variance = rep(0,length(taus))
    for(i in 1:length(taus)){
        tau = taus[i]
        best_variance[i] = get_PVR_big_box_tau(model_name,x,theta,obs,tout,tau)
    }
    plot(taus,best_variance,type='l',xlab='Tau',ylab='Infinite box PVR',cex=1.5,cex.axis=1.5,cex.lab=1.5)
}

plot_PVR_target_goal_tau = function(model_name,x,theta,obs,tout,width=50){
    M = 100
    goals = seq(0,1,length.out=M)[-M]
    results = data.frame(matrix(ncol=3,nrow=0))
    colnames(results) = c("PVR_est","tau","PVR_exact")

    obs_list = matrix(obs,nrow=1)
    lower_limit_list = get_lower(model_name, obs_list,width)
    upper_limit_list = get_upper(model_name, obs_list,width)
    tout_list = c(tout)

    lower_list = get_lower(model_name, obs_list,width)
    upper_list = get_upper(model_name, obs_list,width)

    for(i in 1:length(goals)){
        tau = get_tau(model_name,x,theta,obs,tout,goals[i])
        results[i,] = c(
            get_PVR_big_box_tau(model_name,x,theta,obs,tout,tau),
            tau,
            get_PVR("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,1)$PVR_tout
        )
    }
    plot(goals,results$PVR_exact,type='l',main='Infinite Box PVR',xlab='Target',ylab='Achieved',cex=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,1),ylim=c(0,1))
    lines(goals,results$PVR_est,lty='dashed')
    # abline(a=0,b=1,lty='dashed')
    newresults = results[seq(1,nrow(results),ceiling(M/10)),]
    newgoals = goals[seq(1,nrow(results),ceiling(M/10))]    
    points(newgoals,newresults$PVR_exact)
    text(newgoals,newresults$PVR_exact,paste("tau =",signif(newresults$tau,3)),pos=1,cex=0.8)
    legend(x=0.1,y=0.9,c("Exact","Brownian"),lty = c("solid","dashed"))
}

plot_PVR_target_goal_width = function(model_name,x,theta,obs,tout,tau,width=50){
    M = 100
    goals = seq(0,1,length.out=M)[-M]
    results = data.frame(matrix(ncol=3,nrow=0))
    colnames(results) = c("PVR_est","tau","PVR_exact")

    obs_list = matrix(obs,nrow=1)
    lower_limit_list = get_lower(model_name, obs_list,width)
    upper_limit_list = get_upper(model_name, obs_list,width)
    tout_list = c(tout)

    for(i in 1:length(goals)){
        box = get_box(model_name,theta,tout_list,tau,x,obs_list,goals[i])
        lower_list = box$lower
        upper_list = box$upper
        results[i,] = c(
            get_PVR_big_box_tau(model_name,x,theta,obs,tout,tau),
            tau,
            get_PVR("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,1)$PVR_tout
        )
    }
    plot(goals,results$PVR_exact,type='l',main=paste('Box PVR given tau =',tau),xlab='Target',ylab='Achieved',cex=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,1),ylim=c(0,1))
    lines(c(0,1),c(0,tail(results$PVR_est,1)),lty='dashed')
    # abline(a=0,b=1,lty='dashed')
    newresults = results[seq(1,nrow(results),ceiling(M/10)),]
    newgoals = goals[seq(1,nrow(results),ceiling(M/10))]    
    legend(x=0.1,y=0.9,c("Exact","Brownian"),lty = c("solid","dashed"))
}

plot_PVR_target_goal = function(model_name,x,theta,obs,tout,width=50){
    M = 100
    goals = seq(0,1,length.out=M)[-M]
    results = data.frame(matrix(ncol=3,nrow=0))
    colnames(results) = c("PVR_est","tau","PVR_exact")

    obs_list = matrix(obs,nrow=1)
    lower_limit_list = get_lower(model_name, obs_list,width)
    upper_limit_list = get_upper(model_name, obs_list,width)
    tout_list = c(tout)

    for(i in 1:length(goals)){
        tau = get_tau(model_name,x,theta,obs,tout,sqrt(goals[i]))
        box = get_box(model_name,theta,tout_list,tau,x,obs_list,sqrt(goals[i]))
        lower_list = box$lower
        upper_list = box$upper
        results[i,] = c(
            get_PVR_big_box_tau(model_name,x,theta,obs,tout,tau),
            tau,
            get_PVR("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,1)$PVR_tout
        )
    }
    plot(goals,results$PVR_exact,type='l',main='PVR',xlab='Target',ylab='Achieved',cex=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,1),ylim=c(0,1))
    lines(goals,results$PVR_est*sqrt(goals),lty='dashed')
    # abline(a=0,b=1,lty='dashed')
    newresults = results[seq(1,nrow(results),ceiling(M/10)),]
    newgoals = goals[seq(1,nrow(results),ceiling(M/10))]    
    points(newgoals,newresults$PVR_exact,cex=1)
    text(newgoals,newresults$PVR_exact,paste("tau =",signif(newresults$tau,3)),pos=1,cex=0.8)
    legend(x=0.1,y=0.9,c("Exact","Brownian"),lty = c("solid","dashed"))
}

plot_smooth_box = function(model_name,x,theta,obs,tout,taus,goal,colour){
    l = rep(0,length(taus))
    u = rep(0,length(taus))
    for(i in 1:length(taus)){
            box = get_box_brownian_fast(model_name,theta,tout,taus[i],x,obs,goal)
            l[i] = box$lower
            u[i] = box$upper
    }
    q = c(1,which(diff(l)<0)+1)
    r = c(1,which(diff(u)>0)+1)
    for(i in 1:(length(q))){
            lines(c(tout-taus[q[i]],tout-taus[q[i+1]]),c(l[q[i]],l[q[i+1]]),lwd=2,type='o',col=colour)
    }
    lines(c(0,tout-taus[q[length(q)]]),c(l[q[length(q)]],l[q[length(q)]]),lwd=2,type='o',col=colour)

    for(i in 1:(length(r))){
            lines(c(tout-taus[r[i]],tout-taus[r[i+1]]),c(u[r[i]],u[r[i+1]]),lwd=2,type='o',col=colour)
    }
    lines(c(0,tout-taus[r[length(r)]]),c(u[r[length(r)]],u[r[length(r)]]),lwd=2,type='o',col=colour)
}

plot_smooth_box_changing_goal = function(model_name,x,theta,obs,tout,taus,goal,colour){
    l = rep(0,length(taus))
    u = rep(0,length(taus))
    for(i in 1:length(taus)){
            max_PVR = get_PVR_big_box_tau("BDI",x,theta,obs,tout,tau)
            box = get_box_brownian_fast(model_name,theta,tout,taus[i],x,obs,max(1,goal/max_PVR))
            l[i] = box$lower
            u[i] = box$upper
    }
    q = c(1,which(diff(l)<0)+1)
    r = c(1,which(diff(u)>0)+1)
    for(i in 1:(length(q))){
            lines(c(tout-taus[q[i]],tout-taus[q[i+1]]),c(l[q[i]],l[q[i+1]]),lwd=2,type='o',col=colour)
    }
    lines(c(0,tout-taus[q[length(q)]]),c(l[q[length(q)]],l[q[length(q)]]),lwd=2,type='o',col=colour)

    for(i in 1:(length(r))){
            lines(c(tout-taus[r[i]],tout-taus[r[i+1]]),c(u[r[i]],u[r[i+1]]),lwd=2,type='o',col=colour)
    }
    lines(c(0,tout-taus[r[length(r)]]),c(u[r[length(r)]],u[r[length(r)]]),lwd=2,type='o',col=colour)
}
obs = 110
plot(x=-1,y=-1,xlim=c(0,tout),ylim=c(70,130))
add_obs(tout,obs,'black')
taus = seq(0,tout,length.out=1000)[-1]
plot_smooth_box_changing_goal("BDI",x,theta,obs,tout,taus,0.9,'green')
plot_smooth_box("BDI",x,theta,obs,tout,taus,sqrt(0.9),'black')

get_PVR_big_box_tau("BDI",x,theta,obs,tout,0.4)
0.9/get_PVR_big_box_tau("BDI",x,theta,obs,tout,0.4)
get_box_brownian_fast("BDI",theta,tout,0.4,x,obs,0.9695)
get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,97,116)
get_variance_brownian("BDI",theta,tout,1,x,obs,40,180)

get_q = function(N,M){
    results = rep(0,N)
    for(i in 1:N){
        results[i] = get_likelihood_RB("BDI",x,theta,obs_list,tout_list,M)
    }
    return(results)
}

get_ESJD = function(){

}
