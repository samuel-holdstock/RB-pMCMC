get_chain = function(data,endTime){
    times = data[,1]    
    times = rep(times,each=2)[-1]
    times[length(times)+1] = endTime

    result = matrix(nrow=nrow(data)*2,ncol=ncol(data))
    for(i in 2:ncol(data)){
        result[,i] = rep(data[,i],each=2)
    }
    result[,1] = times
    return(result)
}

plot_chain = function(chain,species,colours,...){
    plot(chain[,c(1,2)],col=colours[1],type='l',...)
    columns = ncol(chain)
    if(columns>2){
        for(i in 3:columns){
            lines(chain[,c(1,i)],col=colours[i-1])
        }
    }
    legend("topleft",legend=species,col=colours,inset=0.05,lty=rep(1,columns-1))
}

add_box = function(lower,upper,tout,tau,colours){
    nboxes = length(lower)
    for(i in 1:nboxes){
        lines(c(tout,tout-tau,tout-tau,tout),c(upper[i],upper[i],lower[i],lower[i]),col=colours[i],lwd=1.5)
    }
}

add_box_goal = function(model_name,goal,theta,x,obs,tout,tau,colours){
    results = get_box_brownian_fast(model_name,theta,tout,tau,x,obs,goal)
    lower = results$lower
    upper = results$upper
    lines(c(tout,tout-tau,tout-tau,tout),c(upper,upper,lower,lower),col=colours,lwd=1.5)
}

add_obs = function(tout,obs,colours){
    npoints = length(obs)
    for(i in 1:npoints){
        points(tout,obs[i],col=colours[i],pch=4,cex=2)
    }
}
