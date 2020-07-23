
discretize <- function(ntinfo_m = NULL, windows = NULL){
  #Round the number of decimals
  ntinfo_m <- round(ntinfo_m, 0)
  
  #Define the interval
  interval <- seq(from = 0, to = 360, by = windows)
  
  #Create list with the intervals
  inter_list <- c()
  
  for(i in 1:length(interval)){
    tryCatch ({
      vect <- (interval[i] + 1):interval[i+1]
      inter_list[[i]] <- vect
    }, error = function(e){
    })
  }
  
  indices <- length(inter_list)
  
  for(j in 1:length(ntinfo_m$alpha)){
    for(i in 1:length(ntinfo_m)){
      la <- which(NA %in% ntinfo_m[j,i])
      if(length(la) == 0){
        for( m in 1:length(inter_list)){
          ind <- which(ntinfo_m[j,i] %in% inter_list[m][[1]])
          if(length(ind) > 0) {
            ntinfo_m[j,i] <- m
          }
        }
      }
    }
  }
  return(ntinfo_m)
}
