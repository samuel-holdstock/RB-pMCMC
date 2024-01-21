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


