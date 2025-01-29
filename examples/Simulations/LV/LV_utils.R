solve_xt_mu = function(model_name,x0,tout,theta,m=10){
  dt = tout/m
  curr_xt = x0
  for(i in 1:m){
    mu = get_mu(model_name,curr_xt,theta)
    curr_xt = curr_xt + mu*dt
  }
  return(curr_xt)
}

get_conditional_xt = function(t,x0,ytout,tout,mu,V,Sigma){
  # A = solve(Sigma*t)
  # B = solve(V+Sigma*(tout-t))
  # C = solve(A+B)
  # normMean = C%*%(A%*%(x0+mu*t)+B%*%(ytout-mu*(tout-t)))
  # normSigma = C
  A = (Sigma*(tout-t))%*%solve(Sigma*tout+V)
  # normMean = x0+mu*(tout-t) + A%*%(ytout-(x0+mu*tout))
  normMean = solve_xt_mu("LV",x0,tout-t,theta,10) + A%*%(ytout-solve_xt_mu("LV",x0,tout,theta,10))
  normSigma = Sigma*(tout-t) - A%*%Sigma*(tout-t)  
  return(list(mean=normMean,covar = normSigma))
}

get_predicted_mean_matrix = function(model_name,tau_list,x,obs_list,tout_list,theta,V=diag(2)*1){
  predicted_mean_matrix = matrix(nrow=length(tout_list),ncol=length(x))
  for(time in 1:length(tout_list)){
    curr_x0 = if(time==1){x}else{obs_list[time-1,]}
    predicted_mean = get_conditional_xt(t=tau_list[time],x0=curr_x0,ytout=obs_list[time,],tout=if(time==1){tout_list[1]}else{tout_list[time]-tout_list[time-1]},mu=get_mu(model_name,(curr_x0+obs_list[time,])/2,theta),V,Sigma=get_covar(model_name,(curr_x0+obs_list[time,])/2,theta))$mean
    predicted_mean_matrix[time,] = predicted_mean
  }
  return(round(predicted_mean_matrix))
}
get_predicted_lower_matrix = function(model_name,tau_list,x,obs_list,tout_list,theta,width,V=diag(2)*1){
  predicted_mean_matrix = round(get_predicted_mean_matrix(model_name,tau_list,x,obs_list,tout_list,theta,V))
  return(get_lower(model_name,predicted_mean_matrix,width))
}

get_predicted_upper_matrix = function(model_name,tau_list,x,obs_list,tout_list,theta,width,V=diag(2)*1){
  predicted_mean_matrix = round(get_predicted_mean_matrix(model_name,tau_list,x,obs_list,tout_list,theta,V))
  return(get_upper(model_name,predicted_mean_matrix,width))
}
plot_box = function(obs_list,time,lower_list,upper_list,col){
  if(ncol(obs_list)==1){
    points(obs_list[1],obs_list[2],col=col,cex=3,lwd=3,pch=16)
  }
  else{
    points(obs_list[time,1],obs_list[time,2],col=col,cex=3,lwd=3,pch=16)
  }
  lines(c(lower_list[time,1],lower_list[time,1]),c(lower_list[time,2],upper_list[time,2]),col=col,lwd=2)
  lines(c(lower_list[time,1],upper_list[time,1]),c(lower_list[time,2],lower_list[time,2]),col=col,lwd=2)
  lines(c(upper_list[time,1],upper_list[time,1]),c(lower_list[time,2],upper_list[time,2]),col=col,lwd=2)
  lines(c(lower_list[time,1],upper_list[time,1]),c(upper_list[time,2],upper_list[time,2]),col=col,lwd=2)
}
box_details = function(d,tout_list){
  print(c("Hit","Stay|Hit"))
  for(i in 1:length(tout_list)){
    print(c(
      mean(d[[i]]$hitS), mean(d[[i]]$inS)/mean(d[[i]]$hitS)
    ))
  }
}
get_box_details = function(d,tout_list){
  M = matrix(ncol=2,nrow=length(d))
  colnames(M) = c("Hit","Stay|Hit")
  for(i in 1:length(tout_list)){
    M[i,] = c(
      mean(d[[i]]$hitS), mean(d[[i]]$inS)/mean(d[[i]]$hitS)
    )
  }
  return(M)
}
get_box_details_string = function(d,tout_list){
  M = list()
  M[[1]] = paste0("Prob Hit\n")
  for(i in 1:length(tout_list)){
    M[[1]] = paste0(M[[1]],"\nObs ",i,": ",
      signif(mean(d[[i]]$hitS),2)
    )
  }
  M[[2]] = paste0("Prob Stay|Hit\n")
  for(i in 1:length(tout_list)){
    M[[2]] = paste0(M[[2]],"\nObs ",i,": ",
      signif(mean(d[[i]]$inS)/mean(d[[i]]$hitS),2)
    )
  }
  return(M)
}

