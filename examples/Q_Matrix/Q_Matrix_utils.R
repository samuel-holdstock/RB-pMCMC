test_exp_time = function(model,lower,upper,theta,S,obs_index,tau,sparse=FALSE){
    Q = get_coffin_matrix(model,lower,upper,theta)
    if(sparse){
        Q = as(Q,"sparseMatrix")
    }
    v = matrix(rep(0,nrow(Q)),nrow=1)
    v[obs_index] = 1
    return(v_exp_Q(v,Q*tau,prec=1e-15,T,T,F))
}
benchmark_exp_time = function(model,lower_mat,upper_mat,theta,S,obs_index,tau,sparse=FALSE){
    n_tests = nrow(lower_mat)
    med_times = rep(0,n_tests)
    Q_size = rep(0,n_tests)
    for(i in 1:n_tests){
        lower = lower_mat[i,]
        upper = upper_mat[i,]
        med_times[i] = summary(microbenchmark(test_exp_time(model,lower,upper,theta,S,obs_index,tau,sparse),unit='micro'))$median
        Q_size[i] = nrow(get_coffin_matrix(model,lower,upper,theta))^2
    }
    return(list(size=Q_size,times=med_times))
}

