### 1D lattice after christen1991 paper
### open boundary conditions and conservative perturbation

### measurement of size s, lifetime t and linear size l for scaling exponent fiiting

### functions:

perform_avalanches <- function(lattice){
  #print("avalanche FUNC:")
  avalanche_lifetime <- 0
  avalanche_size <- 0
  avalance_lin_size <- 0
  r_0 <- which(lattice > z_crit) # initial pos. of avalanche
  #print("R_0")
  #print(r_0)
  
  while (any((lattice > z_crit))) {
    avalanche_lifetime <- avalanche_lifetime + 1
    avalanche_size <- avalanche_size + sum(as.numeric((lattice > z_crit)))
    avalance_lin_size <- max(avalance_lin_size, (which(lattice > z_crit) - r_0))
    #print("actual distance")
    #print(which(lattice > z_crit) - r_0)
    
    tmp_lattice  <- rep(0L, lattice_size)    # auxilliary lattice for calculation of avalanche effects
    
    for(i in 1:lattice_size){ 
      
      if(lattice[i] > z_crit){
        if(i != lattice_size){
          tmp_lattice[i] <- tmp_lattice[i] - 2
          tmp_lattice[i+1] <- tmp_lattice[i+1] + 1
        }
        else {
          tmp_lattice[i] <- tmp_lattice[i] - 1
        }
        
        if(i != 1){
          tmp_lattice[i-1] <- tmp_lattice[i-1] + 1     
        }
        
      }}
    ## add the tmp_lattice to lattice in order to update the lattice
    lattice <- lattice + tmp_lattice
  }
  
  
  
  return (c(lattice, avalanche_size, avalanche_lifetime, avalance_lin_size))
}
  
#### algorithm  
lattice_size <- 20
lattice <- rep(0L, lattice_size)  # main lattice; will contain the z-values
time_steps <- 1e6
z_crit <- 3

s_values <- c()  # list of avalanche sizes
t_values <- c()  # list of avalance lifetimes
l_values <- c()  # list of avalance linear sizes



for(tau in 1:time_steps){
  ## loop until all z<=z_crit  
  ## fill the aux. tmp_lattice conditional on z_crit
  
  tmp <- perform_avalanches(lattice)
  lattice <- tmp[1:lattice_size]
  if(tmp[lattice_size + 1] > 0){   # this cond. checks whether there was an avalanche at all
    s_values[[length(s_values) + 1]] <- tmp[lattice_size + 1]  
    t_values[[length(t_values) + 1]] <- tmp[lattice_size + 2] 
    l_values[[length(l_values) + 1]] <- tmp[lattice_size + 3]  
  }
  
  
    
  ## randomly perturb system via conservative perturbation
  x_rand <- sample(1:lattice_size, 1) 
  
  lattice[x_rand] <- lattice[x_rand] + 1
  if(x_rand > 1){
    lattice[x_rand - 1] <- lattice[x_rand - 1] - 1
  }
  
  
  ## final avalance relaxation
  if(tau == time_steps){
    lattice <- perform_avalanches(lattice)
  }
}


hist(s_values, breaks = 20)
hist(t_values)
hist(l_values)
min(s_values)

sum((s_values - t_values) )

