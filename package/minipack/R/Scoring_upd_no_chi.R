##Newly updated scoring function

scoring_up_wochi <- function(evaluate, bd = bayesian_data){
  
  windows <- bd$windows
  alpha <- bd$alpha
  beta <- bd$beta
  gamma <- bd$gamma
  delta <- bd$delta
  epsilon <- bd$epsilon
  zeta <- bd$zeta
  
  eval <- discretize(ntinfo_m = evaluate, windows = windows)
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  
  amb <- c(alpha[[eval[,"alpha"]]]$beta[eval[,"beta"]] * alpha[[eval[,"alpha"]]]$total / beta[[eval[,"beta"]]]$total,
           alpha[[eval[,"alpha"]]]$gamma[eval[,"gamma"]] * alpha[[eval[,"alpha"]]]$total / gamma[[eval[,"gamma"]]]$total,
           alpha[[eval[,"alpha"]]]$delta[eval[,"delta"]] * alpha[[eval[,"alpha"]]]$total / delta[[eval[,"delta"]]]$total,
           alpha[[eval[,"alpha"]]]$epsilon[eval[,"epsilon"]] * alpha[[eval[,"alpha"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
           alpha[[eval[,"alpha"]]]$zeta[eval[,"zeta"]] * alpha[[eval[,"alpha"]]]$total / zeta[[eval[,"zeta"]]]$total,
           
           
           
           beta[[eval[,"beta"]]]$alpha[eval[,"alpha"]] * beta[[eval[,"beta"]]]$total / alpha[[eval[,"alpha"]]]$total,
           beta[[eval[,"beta"]]]$gamma[eval[,"gamma"]] * beta[[eval[,"beta"]]]$total / gamma[[eval[,"gamma"]]]$total,
           beta[[eval[,"beta"]]]$delta[eval[,"delta"]] * beta[[eval[,"beta"]]]$total / delta[[eval[,"delta"]]]$total,
           beta[[eval[,"beta"]]]$epsilon[eval[,"epsilon"]] * beta[[eval[,"beta"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
           beta[[eval[,"beta"]]]$zeta[eval[,"zeta"]] * beta[[eval[,"beta"]]]$total / zeta[[eval[,"zeta"]]]$total,

           
           
           gamma[[eval[,"gamma"]]]$alpha[eval[,"alpha"]] * gamma[[eval[,"gamma"]]]$total / alpha[[eval[,"alpha"]]]$total,
           gamma[[eval[,"gamma"]]]$beta[eval[,"beta"]] * gamma[[eval[,"gamma"]]]$total / beta[[eval[,"beta"]]]$total,
           gamma[[eval[,"gamma"]]]$delta[eval[,"delta"]] * gamma[[eval[,"gamma"]]]$total / delta[[eval[,"delta"]]]$total,
           gamma[[eval[,"gamma"]]]$epsilon[eval[,"epsilon"]] * gamma[[eval[,"gamma"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
           gamma[[eval[,"gamma"]]]$zeta[eval[,"zeta"]] * gamma[[eval[,"gamma"]]]$total / zeta[[eval[,"zeta"]]]$total,
           
           delta[[eval[,"delta"]]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]]]$total / alpha[[eval[,"alpha"]]]$total,
           delta[[eval[,"delta"]]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]]]$total / beta[[eval[,"beta"]]]$total,
           delta[[eval[,"delta"]]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]]]$total / gamma[[eval[,"gamma"]]]$total,
           delta[[eval[,"delta"]]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
           delta[[eval[,"delta"]]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]]]$total / zeta[[eval[,"zeta"]]]$total,
           
           epsilon[[eval[,"epsilon"]]]$alpha[eval[,"alpha"]] * epsilon[[eval[,"epsilon"]]]$total / alpha[[eval[,"alpha"]]]$total,
           epsilon[[eval[,"epsilon"]]]$beta[eval[,"beta"]] * epsilon[[eval[,"epsilon"]]]$total / beta[[eval[,"beta"]]]$total,
           epsilon[[eval[,"epsilon"]]]$gamma[eval[,"gamma"]] * epsilon[[eval[,"epsilon"]]]$total / gamma[[eval[,"gamma"]]]$total,
           epsilon[[eval[,"epsilon"]]]$delta[eval[,"delta"]] * epsilon[[eval[,"epsilon"]]]$total / delta[[eval[,"delta"]]]$total,
           epsilon[[eval[,"epsilon"]]]$zeta[eval[,"zeta"]] * epsilon[[eval[,"epsilon"]]]$total / zeta[[eval[,"zeta"]]]$total,

           
           zeta[[eval[,"zeta"]]]$alpha[eval[,"alpha"]] * zeta[[eval[,"zeta"]]]$total / alpha[[eval[,"alpha"]]]$total,
           zeta[[eval[,"zeta"]]]$beta[eval[,"beta"]] * zeta[[eval[,"zeta"]]]$total / beta[[eval[,"beta"]]]$total,
           zeta[[eval[,"zeta"]]]$gamma[eval[,"gamma"]] * zeta[[eval[,"zeta"]]]$total / gamma[[eval[,"gamma"]]]$total,
           zeta[[eval[,"zeta"]]]$delta[eval[,"delta"]] * zeta[[eval[,"zeta"]]]$total / delta[[eval[,"delta"]]]$total,
           zeta[[eval[,"zeta"]]]$epsilon[eval[,"epsilon"]] * zeta[[eval[,"zeta"]]]$total / epsilon[[eval[,"epsilon"]]]$total)
  
  ##Process for alpha
  
  #Sum
  
  if(is.na(eval[,"alpha"])){
    amt <- amb
  }else if(eval[,"alpha"]+1 > 36){
    amt <- append(amb, c(alpha[[1]]$beta[eval[,"beta"]] * alpha[[1]]$total / beta[[eval[,"beta"]]]$total,
    alpha[[1]]$gamma[eval[,"gamma"]] * alpha[[1]]$total / gamma[[eval[,"gamma"]]]$total,
    alpha[[1]]$delta[eval[,"delta"]] * alpha[[1]]$total / delta[[eval[,"delta"]]]$total,
    alpha[[1]]$epsilon[eval[,"epsilon"]] * alpha[[1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    alpha[[1]]$zeta[eval[,"zeta"]] * alpha[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"alpha"]+1 <= 36){
    amt <- append(amb, c(alpha[[eval[,"alpha"]+1]]$beta[eval[,"beta"]] * alpha[[eval[,"alpha"]+1]]$total / beta[[eval[,"beta"]]]$total,
    alpha[[eval[,"alpha"]+1]]$gamma[eval[,"gamma"]] * alpha[[eval[,"alpha"]+1]]$total / gamma[[eval[,"gamma"]]]$total,
    alpha[[eval[,"alpha"]+1]]$delta[eval[,"delta"]] * alpha[[eval[,"alpha"]+1]]$total / delta[[eval[,"delta"]]]$total,
    alpha[[eval[,"alpha"]+1]]$epsilon[eval[,"epsilon"]] * alpha[[eval[,"alpha"]+1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    alpha[[eval[,"alpha"]+1]]$zeta[eval[,"zeta"]] * alpha[[eval[,"alpha"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
    
    }
  
  #Substract
  
  if(is.na(eval[,"alpha"])){
    amt <- amb
  }else if(eval[,"alpha"]-1 < 1){
    amt <- append(amt, c(alpha[[36]]$beta[eval[,"beta"]] * alpha[[36]]$total / beta[[eval[,"beta"]]]$total,
    alpha[[36]]$gamma[eval[,"gamma"]] * alpha[[36]]$total / gamma[[eval[,"gamma"]]]$total,
    alpha[[36]]$delta[eval[,"delta"]] * alpha[[36]]$total / delta[[eval[,"delta"]]]$total,
    alpha[[36]]$epsilon[eval[,"epsilon"]] * alpha[[36]]$total / epsilon[[eval[,"epsilon"]]]$total,
    alpha[[36]]$zeta[eval[,"zeta"]] * alpha[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"alpha"]-1 >=1){
    amt <- append(amt, c(alpha[[eval[,"alpha"]-1]]$beta[eval[,"beta"]] * alpha[[eval[,"alpha"]-1]]$total / beta[[eval[,"beta"]]]$total,
    alpha[[eval[,"alpha"]-1]]$gamma[eval[,"gamma"]] * alpha[[eval[,"alpha"]-1]]$total / gamma[[eval[,"gamma"]]]$total,
    alpha[[eval[,"alpha"]-1]]$delta[eval[,"delta"]] * alpha[[eval[,"alpha"]-1]]$total / delta[[eval[,"delta"]]]$total,
    alpha[[eval[,"alpha"]-1]]$epsilon[eval[,"epsilon"]] * alpha[[eval[,"alpha"]-1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    alpha[[eval[,"alpha"]-1]]$zeta[eval[,"zeta"]] * alpha[[eval[,"alpha"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  
  #Process for beta
  
  #Sum
  if(is.na(eval[,"beta"])){
    amt <- amb
  }else if(eval[,"beta"]+1 >36){
    amt <- append(amt,c(beta[[1]]$alpha[eval[,"alpha"]] * beta[[1]]$total / alpha[[eval[,"alpha"]]]$total,
    beta[[1]]$gamma[eval[,"gamma"]] * beta[[1]]$total / gamma[[eval[,"gamma"]]]$total,
    beta[[1]]$delta[eval[,"delta"]] * beta[[1]]$total / delta[[eval[,"delta"]]]$total,
    beta[[1]]$epsilon[eval[,"epsilon"]] * beta[[1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    beta[[1]]$zeta[eval[,"zeta"]] * beta[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"beta"]+1 <= 36){
    amt <- append(amt, c(beta[[eval[,"beta"]+1]]$alpha[eval[,"alpha"]] * beta[[eval[,"beta"]+1]]$total / alpha[[eval[,"alpha"]]]$total,
    beta[[eval[,"beta"]+1]]$gamma[eval[,"gamma"]] * beta[[eval[,"beta"]+1]]$total / gamma[[eval[,"gamma"]]]$total,
    beta[[eval[,"beta"]+1]]$delta[eval[,"delta"]] * beta[[eval[,"beta"]+1]]$total / delta[[eval[,"delta"]]]$total,
    beta[[eval[,"beta"]+1]]$epsilon[eval[,"epsilon"]] * beta[[eval[,"beta"]+1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    beta[[eval[,"beta"]+1]]$zeta[eval[,"zeta"]] * beta[[eval[,"beta"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  
  #Substract
  if(is.na(eval[,"beta"])){
    amt <- amb
  }else if(eval[,"beta"]-1 < 1){
    amt <- append(amt, c(beta[[36]]$alpha[eval[,"alpha"]] * beta[[36]]$total / alpha[[eval[,"alpha"]]]$total,
    beta[[36]]$gamma[eval[,"gamma"]] * beta[[36]]$total / gamma[[eval[,"gamma"]]]$total,
    beta[[36]]$delta[eval[,"delta"]] * beta[[36]]$total / delta[[eval[,"delta"]]]$total,
    beta[[36]]$epsilon[eval[,"epsilon"]] * beta[[36]]$total / epsilon[[eval[,"epsilon"]]]$total,
    beta[[36]]$zeta[eval[,"zeta"]] * beta[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"beta"]-1 >=1){
    amt <- append(amt, c(beta[[eval[,"beta"]-1]]$alpha[eval[,"alpha"]] * beta[[eval[,"beta"]-1]]$total / alpha[[eval[,"alpha"]]]$total,
    beta[[eval[,"beta"]-1]]$gamma[eval[,"gamma"]] * beta[[eval[,"beta"]-1]]$total / gamma[[eval[,"gamma"]]]$total,
    beta[[eval[,"beta"]-1]]$delta[eval[,"delta"]] * beta[[eval[,"beta"]-1]]$total / delta[[eval[,"delta"]]]$total,
    beta[[eval[,"beta"]-1]]$epsilon[eval[,"epsilon"]] * beta[[eval[,"beta"]-1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    beta[[eval[,"beta"]-1]]$zeta[eval[,"zeta"]] * beta[[eval[,"beta"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
    
    
  
  #Process for gamma
  
  #Sum
  
  if(is.na(eval[,"gamma"])){
    amt <- amb
  }else if(eval[,"gamma"] +1 >36){
    amt <- append(amt, c(gamma[[1]]$alpha[eval[,"alpha"]] * gamma[[1]]$total / alpha[[eval[,"alpha"]]]$total,
    gamma[[1]]$beta[eval[,"beta"]] * gamma[[1]]$total / beta[[eval[,"beta"]]]$total,
    gamma[[1]]$delta[eval[,"delta"]] * gamma[[1]]$total / delta[[eval[,"delta"]]]$total,
    gamma[[1]]$epsilon[eval[,"epsilon"]] * gamma[[1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    gamma[[1]]$zeta[eval[,"zeta"]] * gamma[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"gamma"]+1 <= 36){
    amt <- append(amt, c(gamma[[eval[,"gamma"]+1]]$alpha[eval[,"alpha"]] * gamma[[eval[,"gamma"]+1]]$total / alpha[[eval[,"alpha"]]]$total,
    gamma[[eval[,"gamma"]+1]]$beta[eval[,"beta"]] * gamma[[eval[,"gamma"]+1]]$total / beta[[eval[,"beta"]]]$total,
    gamma[[eval[,"gamma"]+1]]$delta[eval[,"delta"]] * gamma[[eval[,"gamma"]+1]]$total / delta[[eval[,"delta"]]]$total,
    gamma[[eval[,"gamma"]+1]]$epsilon[eval[,"epsilon"]] * gamma[[eval[,"gamma"]+1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    gamma[[eval[,"gamma"]+1]]$zeta[eval[,"zeta"]] * gamma[[eval[,"gamma"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  #Substract
  if(is.na(eval[,"gamma"])){
    amt <- amb
  }else if(eval[,"gamma"]-1 < 1){
    amt <- append(amt, c(gamma[[36]]$alpha[eval[,"alpha"]] * gamma[[36]]$total / alpha[[eval[,"alpha"]]]$total,
    gamma[[36]]$beta[eval[,"beta"]] * gamma[[36]]$total / beta[[eval[,"beta"]]]$total,
    gamma[[36]]$delta[eval[,"delta"]] * gamma[[36]]$total / delta[[eval[,"delta"]]]$total,
    gamma[[36]]$epsilon[eval[,"epsilon"]] * gamma[[36]]$total / epsilon[[eval[,"epsilon"]]]$total,
    gamma[[36]]$zeta[eval[,"zeta"]] * gamma[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"gamma"]-1 >=1){
    amt <- append(amt, c(gamma[[eval[,"gamma"]-1]]$alpha[eval[,"alpha"]] * gamma[[eval[,"gamma"]-1]]$total / alpha[[eval[,"alpha"]]]$total,
    gamma[[eval[,"gamma"]-1]]$beta[eval[,"beta"]] * gamma[[eval[,"gamma"]-1]]$total / beta[[eval[,"beta"]]]$total,
    gamma[[eval[,"gamma"]-1]]$delta[eval[,"delta"]] * gamma[[eval[,"gamma"]-1]]$total / delta[[eval[,"delta"]]]$total,
    gamma[[eval[,"gamma"]-1]]$epsilon[eval[,"epsilon"]] * gamma[[eval[,"gamma"]-1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    gamma[[eval[,"gamma"]-1]]$zeta[eval[,"zeta"]] * gamma[[eval[,"gamma"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
    
  
  #Process for delta
  
  #Sum
  if(is.na(eval[,"delta"])){
    amt <- amb
  }else if(eval[,"delta"] +1 >36){
    amt <- append(amt, c(delta[[1]]$alpha[eval[,"alpha"]] * delta[[1]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[1]]$beta[eval[,"beta"]] * delta[[1]]$total / beta[[eval[,"beta"]]]$total,
    delta[[1]]$gamma[eval[,"gamma"]] * delta[[1]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[1]]$epsilon[eval[,"epsilon"]] * delta[[1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[1]]$zeta[eval[,"zeta"]] * delta[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"delta"]+1 <= 36){
    amt <- append(amt, c(delta[[eval[,"delta"]+1]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]+1]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[eval[,"delta"]+1]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]+1]]$total / beta[[eval[,"beta"]]]$total,
    delta[[eval[,"delta"]+1]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]+1]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[eval[,"delta"]+1]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]+1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[eval[,"delta"]+1]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  #Substract
  
  if(is.na(eval[,"delta"])){
    amt <- amb
  }else if(eval[,"delta"]-1 < 1){
    amt <- append(amt, c(delta[[36]]$alpha[eval[,"alpha"]] * delta[[36]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[36]]$beta[eval[,"beta"]] * delta[[36]]$total / beta[[eval[,"beta"]]]$total,
    delta[[36]]$gamma[eval[,"gamma"]] * delta[[36]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[36]]$epsilon[eval[,"epsilon"]] * delta[[36]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[36]]$zeta[eval[,"zeta"]] * delta[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"delta"]-1 >=1){
    amt <- append(amt,  c(delta[[eval[,"delta"]-1]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]-1]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[eval[,"delta"]-1]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]-1]]$total / beta[[eval[,"beta"]]]$total,
    delta[[eval[,"delta"]-1]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]-1]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[eval[,"delta"]-1]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]-1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[eval[,"delta"]-1]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
    
    
  
  
  
  #Process for epsilon
  
  #Sum
  if(is.na(eval[,"epsilon"])){
    amt <- amb
  }else if(eval[,"epsilon"] +1 >36){
    amt <- append(amt, c(epsilon[[1]]$alpha[eval[,"alpha"]] * epsilon[[1]]$total / alpha[[eval[,"alpha"]]]$total,
    epsilon[[1]]$beta[eval[,"beta"]] * epsilon[[1]]$total / beta[[eval[,"beta"]]]$total,
    epsilon[[1]]$gamma[eval[,"gamma"]] * epsilon[[1]]$total / gamma[[eval[,"gamma"]]]$total,
    epsilon[[1]]$delta[eval[,"delta"]] * epsilon[[1]]$total / delta[[eval[,"delta"]]]$total,
    epsilon[[1]]$zeta[eval[,"zeta"]] * epsilon[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"epsilon"]+1 <= 36){
    amt <- append(amt, c(epsilon[[eval[,"epsilon"]+1]]$alpha[eval[,"alpha"]] * epsilon[[eval[,"epsilon"]+1]]$total / alpha[[eval[,"alpha"]]]$total,
    epsilon[[eval[,"epsilon"]+1]]$beta[eval[,"beta"]] * epsilon[[eval[,"epsilon"]+1]]$total / beta[[eval[,"beta"]]]$total,
    epsilon[[eval[,"epsilon"]+1]]$gamma[eval[,"gamma"]] * epsilon[[eval[,"epsilon"]+1]]$total / gamma[[eval[,"gamma"]]]$total,
    epsilon[[eval[,"epsilon"]+1]]$delta[eval[,"delta"]] * epsilon[[eval[,"epsilon"]+1]]$total / delta[[eval[,"delta"]]]$total,
    epsilon[[eval[,"epsilon"]+1]]$zeta[eval[,"zeta"]] * epsilon[[eval[,"epsilon"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  #Substract
  if(is.na(eval[,"epsilon"])){
    amt <- amb
  }else if(eval[,"epsilon"]-1 < 1){
    amt <- append(amt, c(epsilon[[36]]$alpha[eval[,"alpha"]] * epsilon[[36]]$total / alpha[[eval[,"alpha"]]]$total,
    epsilon[[36]]$beta[eval[,"beta"]] * epsilon[[36]]$total / beta[[eval[,"beta"]]]$total,
    epsilon[[36]]$gamma[eval[,"gamma"]] * epsilon[[36]]$total / gamma[[eval[,"gamma"]]]$total,
    epsilon[[36]]$delta[eval[,"delta"]] * epsilon[[36]]$total / delta[[eval[,"delta"]]]$total,
    epsilon[[36]]$zeta[eval[,"zeta"]] * epsilon[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"epsilon"]-1 >=1){
    amt <- append(amt, c(epsilon[[eval[,"epsilon"]-1]]$alpha[eval[,"alpha"]] * epsilon[[eval[,"epsilon"]-1]]$total / alpha[[eval[,"alpha"]]]$total,
    epsilon[[eval[,"epsilon"]-1]]$beta[eval[,"beta"]] * epsilon[[eval[,"epsilon"]-1]]$total / beta[[eval[,"beta"]]]$total,
    epsilon[[eval[,"epsilon"]-1]]$gamma[eval[,"gamma"]] * epsilon[[eval[,"epsilon"]-1]]$total / gamma[[eval[,"gamma"]]]$total,
    epsilon[[eval[,"epsilon"]-1]]$delta[eval[,"delta"]] * epsilon[[eval[,"epsilon"]-1]]$total / delta[[eval[,"delta"]]]$total,
    epsilon[[eval[,"epsilon"]-1]]$zeta[eval[,"zeta"]] * epsilon[[eval[,"epsilon"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  
  #Process for zeta
  
  #Sum
  if(is.na(eval[,"zeta"])){
    amt <- amb
  }else if(eval[,"zeta"] +1 >36){
    amt <- append(amt, c(zeta[[1]]$alpha[eval[,"alpha"]] * zeta[[1]]$total / alpha[[eval[,"alpha"]]]$total,
                         zeta[[1]]$beta[eval[,"beta"]] * zeta[[1]]$total / beta[[eval[,"beta"]]]$total,
                         zeta[[1]]$gamma[eval[,"gamma"]] * zeta[[1]]$total / gamma[[eval[,"gamma"]]]$total,
                         zeta[[1]]$delta[eval[,"delta"]] * zeta[[1]]$total / delta[[eval[,"delta"]]]$total,
                         zeta[[1]]$zeta[eval[,"zeta"]] * zeta[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"zeta"]+1 <= 36){
    amt <- append(amt, c(zeta[[eval[,"zeta"]+1]]$alpha[eval[,"alpha"]] * zeta[[eval[,"zeta"]+1]]$total / alpha[[eval[,"alpha"]]]$total,
                         zeta[[eval[,"zeta"]+1]]$beta[eval[,"beta"]] * zeta[[eval[,"zeta"]+1]]$total / beta[[eval[,"beta"]]]$total,
                         zeta[[eval[,"zeta"]+1]]$gamma[eval[,"gamma"]] * zeta[[eval[,"zeta"]+1]]$total / gamma[[eval[,"gamma"]]]$total,
                         zeta[[eval[,"zeta"]+1]]$delta[eval[,"delta"]] * zeta[[eval[,"zeta"]+1]]$total / delta[[eval[,"delta"]]]$total,
                         zeta[[eval[,"zeta"]+1]]$zeta[eval[,"zeta"]] * zeta[[eval[,"zeta"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  #Substract
  
  if(is.na(eval[,"zeta"])){
    amt <- amb
  }else if(eval[,"zeta"]-1 < 1){
    amt <- append(amt, c(zeta[[36]]$alpha[eval[,"alpha"]] * zeta[[36]]$total / alpha[[eval[,"alpha"]]]$total,
                         zeta[[36]]$beta[eval[,"beta"]] * zeta[[36]]$total / beta[[eval[,"beta"]]]$total,
                         zeta[[36]]$gamma[eval[,"gamma"]] * zeta[[36]]$total / gamma[[eval[,"gamma"]]]$total,
                         zeta[[36]]$delta[eval[,"delta"]] * epsilon[[36]]$total / delta[[eval[,"delta"]]]$total,
                         zeta[[36]]$zeta[eval[,"zeta"]] * zeta[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"zeta"]-1 >=1){
    amt <- append(amt, c(zeta[[eval[,"zeta"]-1]]$alpha[eval[,"alpha"]] * zeta[[eval[,"zeta"]-1]]$total / alpha[[eval[,"alpha"]]]$total,
                         zeta[[eval[,"zeta"]-1]]$beta[eval[,"beta"]] * zeta[[eval[,"zeta"]-1]]$total / beta[[eval[,"beta"]]]$total,
                         zeta[[eval[,"zeta"]-1]]$gamma[eval[,"gamma"]] * zeta[[eval[,"zeta"]-1]]$total / gamma[[eval[,"gamma"]]]$total,
                         zeta[[eval[,"zeta"]-1]]$delta[eval[,"delta"]] * zeta[[eval[,"zeta"]-1]]$total / delta[[eval[,"delta"]]]$total,
                         zeta[[eval[,"zeta"]-1]]$zeta[eval[,"zeta"]] * zeta[[eval[,"zeta"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  return(round(sum(amt)/length(amt), 3))
}
