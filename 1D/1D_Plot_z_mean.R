### 1D lattice after christen1991 paper
### open boundary conditions and conservative perturbation


### functions:

perform_avalanches <- function(lattice){
  while (any((lattice > z_crit))) {
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
  return (lattice)
}
  
#### algorithm  
lattice_size <- 20

time_steps <- 3500
z_crit_max <- 6
color <- c("green4","blue4","red2","orange2","pink2")


for(z_crit in z_crit_max:1){
  lattice <- rep(0L, lattice_size)  # main lattice; will contain the z-values
  z_means <- numeric(time_steps)
  
  for(t in 1:time_steps){
    
    lattice <- perform_avalanches(lattice)
    
    
    ## randomly perturb system via conservative perturbation
    x_rand <- sample(1:lattice_size, 1) 
    
    lattice[x_rand] <- lattice[x_rand] + 1
    if(x_rand > 1){
      lattice[x_rand - 1] <- lattice[x_rand - 1] - 1
    }
    
    ## compute z_mean[t]
    z_means[t] <- mean(lattice)
    
    
    ## final avalance relaxation
    if(t == time_steps){
      lattice <- perform_avalanches(lattice)
    }
  }
  if(z_crit == z_crit_max){
    plot(z_means,type = "l", xlab = expression(paste("time steps ", tau)), 
         ylab = expression(paste("<z(", tau, ")>")))  
  } else {
    lines(z_means)
    lines(z_means, col=color[z_crit])
    
  }
  
}
  
legend(0, 6, legend=c(expression(paste('z'['crit'], ' = 6')),
                      expression(paste('z'['crit'], ' = 5')),
                      expression(paste('z'['crit'], ' = 4')),
                      expression(paste('z'['crit'], ' = 3')),
                      expression(paste('z'['crit'], ' = 2')),
                      expression(paste('z'['crit'], ' = 1'))),
       col=c("black","pink2","orange2","red2","blue4","green4"), lty=1:1, cex=0.75)



## PLOTS
# barplot(lattice)
# plot(z_means, type = "l")
