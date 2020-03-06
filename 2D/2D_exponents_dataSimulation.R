### 2D lattice after christen1991 paper
### open boundary conditions and conservative perturbation

### functions:
perform_avalanches <- function(lattice){
  #print(" ")
  #print("avalanche func start:")
  
  avalanche_lifetime <- 0
  avalanche_size <- 0
  avalance_lin_size <- 0
  
  
  r_0 <- which(lattice > z_crit, arr.ind = TRUE) # initial pos. of avalanche 
  #print(c("r_0:", r_0))
  
  while (any((lattice > z_crit))) {
    avalanche_lifetime <- avalanche_lifetime + 1
    avalanche_size <- avalanche_size + sum(as.numeric((lattice > z_crit)))
    
    
    
    tmp_lattice  <- matrix(0L, nrow = lattice_size, ncol = lattice_size)    # auxilliary lattice for calculation of avalanche effects
    
    #which(lattice > z_crit) %% lattice_size   # row number
    #(which(lattice > z_crit) %/% lattice_size) + 1 # column number
    
    for(i in 1:lattice_size){ 
      for(j in 1:lattice_size){  
        if(lattice[i,j] > z_crit){
          
          actual_lin_dis <- sqrt((i - r_0[1])^2 + (j - r_0[2])^2)# distance between r_0 actual critical point
          
          # 1 is min. lin size; 1.21 ist min. size wenn s = 2 
          if(avalanche_size > 1 & actual_lin_dis < 1.21) {actual_lin_dis = 1.21}

          avalance_lin_size <- max(1, max(avalance_lin_size, actual_lin_dis))
          

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
  
    
  
  }
  #print(c("s:", avalanche_size, "t:", avalanche_lifetime, "l:", avalance_lin_size))
  return (c(lattice, avalanche_size, avalanche_lifetime, avalance_lin_size))
}


#### algorithm  
library(plot.matrix)
library(randtoolbox)

lattice_size <- 45
z_crit <- 8      # avalanche condtion
time_steps <- 1e7  # 3*1e7 lasts > 25min with lattice_size=10

lattice <- matrix(0, nrow = lattice_size, ncol = lattice_size)     # the main lattice for the simulation
random_gen_choice <- 1 # 1= std. sample(); 2= random package real random numbers


perturbat_conser <- TRUE  # perturbation conservative (TRUE) or non-conservative (FALSE)


s_values <- c()  # list of avalanche sizes
t_values <- c()  # list of avalance lifetimes
l_values <- c()  # list of avalance linear sizes


for(t in 1:time_steps){
  tmp <- perform_avalanches(lattice)
  lattice <- matrix(tmp[1:lattice_size^2], nrow = lattice_size, ncol = lattice_size)
  
  if(tmp[lattice_size^2 + 1] > 0){   # this cond. checks whether there was an avalanche at all
    s_values[[length(s_values) + 1]] <- tmp[lattice_size^2 + 1]  
    t_values[[length(t_values) + 1]] <- tmp[lattice_size^2 + 2] 
    l_values[[length(l_values) + 1]] <- tmp[lattice_size^2 + 3]  
  }
  
  #### randomly place one grain of sand by z<-z+1 on random pos.
  
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
  
  
  if(t %% 1e4 == 0){
    print(c(t/time_steps * 100, "% steps done", 
            "<z>", mean(lattice)))
  }
  
  
  #if(t %% 1000 == 0){
  #  plot(lattice,  breaks=c(0, 1, 2,3,4,5,6,7,8,9), digits = 0)
  #}
  
  
  ## final avalance relaxation
  if(t == time_steps){
    lattice <- perform_avalanches(lattice)
  }
  
  
}

# final step: creationof data frame containing all avalanche events
simulation_data <- data.frame(s=s_values,t=t_values,l=l_values)


