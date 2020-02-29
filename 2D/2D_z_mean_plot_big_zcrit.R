### 2D lattice after christen1991 paper
### open boundary conditions and conservative perturbation

### functions:

perform_avalanches <- function(lattice){
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
lattice_size <- 15
time_steps <- 0.5*1e5  # 3*1e7 lasts > 25min with lattice_size=10
plotting <- F


perturbat_conser <- FALSE  # perturbation conservative (TRUE) or non-conservative (FALSE)

color <- c("green4","blue4","red2")

z_crit_max <- 100
i <- 1

for(z_crit in c(z_crit_max, 75, 50, 10)){
  
  lattice <- matrix(0L, nrow = lattice_size, ncol = lattice_size)     # the main lattice for the simulation
  z_means <- numeric(time_steps)
  
  for(t in 1:time_steps){
    ## loop until all z<=z_crit  
  
      
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
    lines(z_means, col=color[i])
    i <- i + 1
  }
}

grid()
legend(0, 100, legend=c(expression(paste('z'['crit'], ' = 100')),
                      expression(paste('z'['crit'], ' = 75')),
                      expression(paste('z'['crit'], ' = 50')),
                      expression(paste('z'['crit'], ' = 10'))),
       col=c("black","green4","blue4","red2"), lty=1:1, cex=0.8)

