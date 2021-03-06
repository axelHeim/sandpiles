### 2D lattice after christen1991 paper
### open boundary conditions and non-cons. OR conservative perturbation

### functions:

perform_avalanches <- function(lattice){
  if(plotting == TRUE){
    plot(lattice,  breaks=c(0, 1, 2,3,4,5,6), digits = 0)
    Sys.sleep(1)
    
  }
  
  
  while (any((lattice > z_crit))) {
    
    tmp_lattice  <- matrix(0L, nrow = lattice_size, ncol = lattice_size)    # auxilliary lattice for calculation of avalanche effects
    
    #which(lattice > z_crit) %% lattice_size   # row number
    #(which(lattice > z_crit) %/% lattice_size) + 1 # column number
    
    for(i in 1:lattice_size){ 
      for(j in 1:lattice_size){  
        if(lattice[i,j] > z_crit){
          
          tmp_lattice[i,j] <- tmp_lattice[i,j] - 4
          
          if(i != lattice_size){
            tmp_lattice[i+1,j] <- tmp_lattice[i+1,j] + 1
          }
          else {
            tmp_lattice[i,j] <- tmp_lattice[i,j] + 1
          }
          if(i != 1){
            tmp_lattice[i-1,j] <- tmp_lattice[i-1,j] + 1     
          }
          
          
          
          if(j != lattice_size){
            tmp_lattice[i,j+1] <- tmp_lattice[i,j+1] + 1
          }
          else {
            tmp_lattice[i,j] <- tmp_lattice[i,j] + 1
          }
          if(j != 1){
            tmp_lattice[i,j-1] <- tmp_lattice[i,j-1] + 1     
          }
          
          
      }}}
    ## add the tmp_lattice to lattice in order to update the lattice
    lattice <- lattice + tmp_lattice
  
    if(plotting == TRUE){
      plot(lattice,  breaks=c(0, 1, 2,3,4,5,6), digits = 0)
      Sys.sleep(1)
      
    }
    
  
    }
  return (lattice)
}


#### algorithm  
library(plot.matrix)
lattice_size <- 9
z_crit <- 8      # avalanche condtion
time_steps <- 1e6  # 3*1e7 lasts > 25min with lattice_size=10
plotting <- F


perturbat_conser <- FALSE  # perturbation conservative (TRUE) or non-conservative (FALSE)

lattice <- matrix(0L, nrow = lattice_size, ncol = lattice_size)     # the main lattice for the simulation


for(t in 1:time_steps){
  ## loop until all z<=z_crit  
  ## fill the aux. tmp_lattice conditional on z_crit
  #print("timestep")
  #print(t)
    
  lattice <- perform_avalanches(lattice)
  
  
  ## randomly place one grain of sand by z<-z+1 on random pos.
  x_rand <- sample(1:lattice_size, 1) 
  y_rand <- sample(1:lattice_size, 1) 
  
  if(perturbat_conser == TRUE){ # conservative perturbation
    lattice[x_rand, y_rand] <- lattice[x_rand, y_rand] + 2
    if(x_rand > 1){
      lattice[x_rand - 1, y_rand] <- lattice[x_rand - 1, y_rand] - 1
    }
    if(y_rand > 1){
      lattice[x_rand, y_rand - 1] <- lattice[x_rand, y_rand - 1] - 1
    }
  } else { # non-conservative perturbation
    lattice[x_rand, y_rand] <- lattice[x_rand, y_rand] + 1
  }
  
  
  ## final avalance relaxation
  if(t == time_steps){
    lattice <- perform_avalanches(lattice)
  }
  
  if(t %% 1e5 == 0){
    print("% steps done:")
    print(t/time_steps * 100)
    print(c("<z>",mean(lattice)))
  }
}

plot(lattice,  breaks=c(0, 1, 2,3,4,5,6,7,8,9), digits = 0)
sum(as.numeric(lattice == z_crit))/(lattice_size**2)
mean(lattice)

