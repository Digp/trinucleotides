
probability <- function(angle = NULL, ntinfo_disc = NULL){
  
  ntinfo_m <- ntinfo_disc
  len <- max(ntinfo_m, na.rm=TRUE)
  inter_list <- 1:len
  
  ######################Run_Code###################
  
  #Calculate the number of nucleotides of our dataset
  tot <- length(ntinfo_disc$alpha)
  
  ######For alpha
  
  if(angle == "alpha"){
    alpha <- vector("list", length = length(inter_list))
      
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"alpha"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        alpha[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"beta"] == j)
          if (length(indim) == 0) {
            alpha[[i]]$beta[j] <- 0
          } else {
            alpha[[i]]$beta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"gamma"] == j)
          if (length(indim) == 0) {
            alpha[[i]]$gamma[j] <- 0
          } else {
            alpha[[i]]$gamma[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"delta"] == j)
          if (length(indim) == 0) {
            alpha[[i]]$delta[j] <- 0
          } else {
            alpha[[i]]$delta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"epsilon"] == j)
          if (length(indim) == 0) {
            alpha[[i]]$epsilon[j] <- 0
          } else {
            alpha[[i]]$epsilon[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"zeta"] == j)
          if (length(indim) == 0) {
            alpha[[i]]$zeta[j] <- 0
          } else {
            alpha[[i]]$zeta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"chi"] == j)
          if (length(indim) == 0) {
            alpha[[i]]$chi[j] <- 0
          } else {
            alpha[[i]]$chi[j] <- length(indim) / n
          }
        }
      }else{
        alpha[[i]]$beta <- rep(0, len)
        alpha[[i]]$gamma <- rep(0, len)
        alpha[[i]]$delta <- rep(0, len)
        alpha[[i]]$epsilon <- rep(0, len)
        alpha[[i]]$zeta <- rep(0, len)
        alpha[[i]]$chi <- rep(0, len)
        alpha[[i]]$total <- 0
        
      }
    }
    return(alpha)
  }
  
  
  #####For beta
    
  if(angle == "beta"){
    beta <- vector("list", length = length(inter_list))
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"beta"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        beta[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"alpha"] == j)
          if (length(indim) == 0) {
            beta[[i]]$alpha[j] <- 0
          } else {
            beta[[i]]$alpha[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"gamma"] == j)
          if (length(indim) == 0) {
            beta[[i]]$gamma[j] <- 0
          } else {
            beta[[i]]$gamma[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"delta"] == j)
          if (length(indim) == 0) {
            beta[[i]]$delta[j] <- 0
          } else {
            beta[[i]]$delta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"epsilon"] == j)
          if (length(indim) == 0) {
            beta[[i]]$epsilon[j] <- 0
          } else {
            beta[[i]]$epsilon[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"zeta"] == j)
          if (length(indim) == 0) {
            beta[[i]]$zeta[j] <- 0
          } else {
            beta[[i]]$zeta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"chi"] == j)
          if (length(indim) == 0) {
            beta[[i]]$chi[j] <- 0
          } else {
            beta[[i]]$chi[j] <- length(indim) / n
          }
        }
      }else{
        beta[[i]]$alpha <- rep(0, len)
        beta[[i]]$gamma <- rep(0, len)
        beta[[i]]$delta <- rep(0, len)
        beta[[i]]$epsilon <- rep(0, len)
        beta[[i]]$zeta <- rep(0, len)
        beta[[i]]$chi <- rep(0, len)
        beta[[i]]$total <- 0
      }
    }
    return(beta)
  }
  

  
  #####For gamma
  
  if(angle == "gamma"){
    gamma <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"gamma"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        gamma[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"alpha"] == j)
          if (length(indim) == 0) {
            gamma[[i]]$alpha[j] <- 0
          } else {
            gamma[[i]]$alpha[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"beta"] == j)
          if (length(indim) == 0) {
            gamma[[i]]$beta[j] <- 0
          } else {
            gamma[[i]]$beta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"delta"] == j)
          if (length(indim) == 0) {
            gamma[[i]]$delta[j] <- 0
          } else {
            gamma[[i]]$delta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"epsilon"] == j)
          if (length(indim) == 0) {
            gamma[[i]]$epsilon[j] <- 0
          } else {
            gamma[[i]]$epsilon[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"zeta"] == j)
          if (length(indim) == 0) {
            gamma[[i]]$zeta[j] <- 0
          } else {
            gamma[[i]]$zeta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"chi"] == j)
          if (length(indim) == 0) {
            gamma[[i]]$chi[j] <- 0
          } else {
            gamma[[i]]$chi[j] <- length(indim) / n
          }
        }
      }else{
        gamma[[i]]$alpha <- rep(0, len)
        gamma[[i]]$beta <- rep(0, len)
        gamma[[i]]$delta <- rep(0, len)
        gamma[[i]]$epsilon <- rep(0, len)
        gamma[[i]]$zeta <- rep(0, len)
        gamma[[i]]$chi <- rep(0, len)
        gamma[[i]]$total <- 0
        
      }
    }
    return(gamma)
  }
  
  #####For delta
  if(angle == "delta"){
    delta <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"delta"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        delta[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"alpha"] == j)
          if (length(indim) == 0) {
            delta[[i]]$alpha[j] <- 0
          } else {
            delta[[i]]$alpha[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"beta"] == j)
          if (length(indim) == 0) {
            delta[[i]]$beta[j] <- 0
          } else {
            delta[[i]]$beta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"gamma"] == j)
          if (length(indim) == 0) {
            delta[[i]]$gamma[j] <- 0
          } else {
            delta[[i]]$gamma[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"epsilon"] == j)
          if (length(indim) == 0) {
            delta[[i]]$epsilon[j] <- 0
          } else {
            delta[[i]]$epsilon[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"zeta"] == j)
          if (length(indim) == 0) {
            delta[[i]]$zeta[j] <- 0
          } else {
            delta[[i]]$zeta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"chi"] == j)
          if (length(indim) == 0) {
            delta[[i]]$chi[j] <- 0
          } else {
            delta[[i]]$chi[j] <- length(indim) / n
          }
        }
      }else{
        delta[[i]]$alpha <- rep(0, len)
        delta[[i]]$beta <- rep(0, len)
        delta[[i]]$gamma <- rep(0, len)
        delta[[i]]$epsilon <- rep(0, len)
        delta[[i]]$zeta <- rep(0, len)
        delta[[i]]$chi <- rep(0, len)
        delta[[i]]$total <- 0

      }
    }
    return(delta)
  }
  
  ##### For epsilon
  if(angle == "epsilon"){
    epsilon <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"epsilon"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        epsilon[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"alpha"] == j)
          if (length(indim) == 0) {
            epsilon[[i]]$alpha[j] <- 0
          } else {
            epsilon[[i]]$alpha[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"beta"] == j)
          if (length(indim) == 0) {
            epsilon[[i]]$beta[j] <- 0
          } else {
            epsilon[[i]]$beta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"gamma"] == j)
          if (length(indim) == 0) {
            epsilon[[i]]$gamma[j] <- 0
          } else {
            epsilon[[i]]$gamma[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"delta"] == j)
          if (length(indim) == 0) {
            epsilon[[i]]$delta[j] <- 0
          } else {
            epsilon[[i]]$delta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"zeta"] == j)
          if (length(indim) == 0) {
            epsilon[[i]]$zeta[j] <- 0
          } else {
            epsilon[[i]]$zeta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"chi"] == j)
          if (length(indim) == 0) {
            epsilon[[i]]$chi[j] <- 0
          } else {
            epsilon[[i]]$chi[j] <- length(indim) / n
          }
        }
      }else{
        epsilon[[i]]$alpha <- rep(0, len) 
        epsilon[[i]]$beta <- rep(0, len)
        epsilon[[i]]$gamma <- rep(0, len)
        epsilon[[i]]$delta <- rep(0, len)
        epsilon[[i]]$zeta <- rep(0, len)
        epsilon[[i]]$chi <- rep(0, len)
        epsilon[[i]]$total <- 0
        
      }
    }
    return(epsilon)
  }
  
  #####For zeta
  if(angle == "zeta"){
    zeta <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"zeta"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        zeta[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"alpha"] == j)
          if (length(indim) == 0) {
            zeta[[i]]$alpha[j] <- 0
          } else {
            zeta[[i]]$alpha[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"beta"] == j)
          if (length(indim) == 0) {
            zeta[[i]]$beta[j] <- 0
          } else {
            zeta[[i]]$beta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"gamma"] == j)
          if (length(indim) == 0) {
            zeta[[i]]$gamma[j] <- 0
          } else {
            zeta[[i]]$gamma[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"delta"] == j)
          if (length(indim) == 0) {
            zeta[[i]]$delta[j] <- 0
          } else {
            zeta[[i]]$delta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"epsilon"] == j)
          if (length(indim) == 0) {
            zeta[[i]]$epsilon[j] <- 0
          } else {
            zeta[[i]]$epsilon[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"chi"] == j)
          if (length(indim) == 0) {
            zeta[[i]]$chi[j] <- 0
          } else {
            zeta[[i]]$chi[j] <- length(indim) / n
          }
        }
      }else{
        zeta[[i]]$alpha <- rep(0, len)
        zeta[[i]]$beta <- rep(0, len)
        zeta[[i]]$gamma <- rep(0, len)
        zeta[[i]]$delta <- rep(0, len)
        zeta[[i]]$epsilon <- rep(0, len)
        zeta[[i]]$chi <- rep(0, len)
        zeta[[i]]$total <- 0
      }  
    }
    return(zeta)
  }
  
  ####For chi
  if(angle == "chi"){
    chi <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"chi"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        chi[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"alpha"] == j)
          if (length(indim) == 0) {
            chi[[i]]$alpha[j] <- 0
          } else {
            chi[[i]]$alpha[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"beta"] == j)
          if (length(indim) == 0) {
            chi[[i]]$beta[j] <- 0
          } else {
            chi[[i]]$beta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"gamma"] == j)
          if (length(indim) == 0) {
            chi[[i]]$gamma[j] <- 0
          } else {
            chi[[i]]$gamma[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"delta"] == j)
          if (length(indim) == 0) {
            chi[[i]]$delta[j] <- 0
          } else {
            chi[[i]]$delta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"epsilon"] == j)
          if (length(indim) == 0) {
            chi[[i]]$epsilon[j] <- 0
          } else {
            chi[[i]]$epsilon[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"zeta"] == j)
          if (length(indim) == 0) {
            chi[[i]]$zeta[j] <- 0
          } else {
            chi[[i]]$zeta[j] <- length(indim) / n
          }
        }
      }else{
        chi[[i]]$alpha <- rep(0, len)
        chi[[i]]$beta <- rep(0, len)
        chi[[i]]$gamma <- rep(0, len)
        chi[[i]]$delta <- rep(0, len)
        chi[[i]]$epsilon <- rep(0, len)
        chi[[i]]$zeta <- rep(0, len)
        chi[[i]]$total <- 0
      }
    }
    return(chi)
  }
}
