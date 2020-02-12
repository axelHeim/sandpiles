library(plot.matrix)
lattice_size <- 12
z_crit <- 6       # avalanche condtion
time_steps <- 10000


lattice <- matrix(0L, nrow = lattice_size, ncol = lattice_size)     # the main lattice for the simulation


for(t in 1:time_steps){
  ## loop until all z<=z_crit  
  ## fill the aux. tmp_lattice conditional on z_crit
  print("timestep")
  print(t)
    
  while (any((lattice > z_crit))) {
    tmp_lattice  <- matrix(0L, nrow = lattice_size, ncol = lattice_size)    # auxilliary lattice for calculation of avalanches
    
    for(i in 2:(lattice_size - 1)){
      for(j in 2:(lattice_size - 1)){
        if(lattice[i,j] > z_crit){
          tmp_lattice[i,j] <- tmp_lattice[i,j] - 4   
          tmp_lattice[i+1,j] <- tmp_lattice[i+1,j] + 1   
          tmp_lattice[i,j-1] <- tmp_lattice[i,j-1] + 1   
          tmp_lattice[i-1,j] <- tmp_lattice[i-1,j] + 1   
          tmp_lattice[i,j+1] <- tmp_lattice[i,j+1] + 1
        }}}
    ## add the tmp_lattice to lattice in order to update the lattice
    lattice <- lattice + tmp_lattice
    
    ## set boundaries to zero
    for(i in 1:lattice_size){
      lattice[i,1] <- 0
      lattice[i,lattice_size] <- 0
    }
    for(j in 1:lattice_size){
      lattice[1,j] <- 0
      lattice[lattice_size,j] <- 0
    }  
    }
    
  ## randomly place one grain of sand by z<-z+1 on random pos.
  x_rand <- sample(2:(lattice_size - 1), 1) 
  y_rand <- sample(2:(lattice_size - 1), 1) 
  
  lattice[x_rand, y_rand] <- lattice[x_rand, y_rand] + 1
}

plot(lattice)
