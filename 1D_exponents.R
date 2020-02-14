### 1D lattice after christen1991 paper
### open boundary conditions and conservative perturbation

### measurement of size s, lifetime t and linear size l for scaling exponent fiiting

### functions:

perform_avalanches <- function(lattice){
  while (any((lattice > z_crit))) {
    tmp_lattice  <- rep(0L, lattice_size)    # auxilliary lattice for calculation of avalanche effects
    
    for(i in 1:lattice_size){ # skip first entry because of bound. cond.
      
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
  return (lattice)
}
  
#### algorithm  
lattice_size <- 80
lattice <- rep(0L, lattice_size)  # main lattice; will contain the z-values
time_steps <- 10000
z_crit <- 1

  
  
  


for(t in 1:time_steps){
  ## loop until all z<=z_crit  
  ## fill the aux. tmp_lattice conditional on z_crit
  
  lattice <- perform_avalanches(lattice)
    
  ## randomly perturb system via conservative perturbation
  x_rand <- sample(1:lattice_size, 1) 
  
  lattice[x_rand] <- lattice[x_rand] + 1
  if(x_rand > 1){
    lattice[x_rand - 1] <- lattice[x_rand - 1] - 1
  }
  
  
  
  ## final avalance relaxation
  if(t == time_steps){
    lattice <- perform_avalanches(lattice)
  }
}


