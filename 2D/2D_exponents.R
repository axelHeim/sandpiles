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
          avalance_lin_size <- max(avalance_lin_size, actual_lin_dis)
          

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
lattice_size <- 40
z_crit <- 8      # avalanche condtion
time_steps <- 5*1e6  # 3*1e7 lasts > 25min with lattice_size=10

lattice <- matrix(0, nrow = lattice_size, ncol = lattice_size)     # the main lattice for the simulation

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
  
  
  ## randomly place one grain of sand by z<-z+1 on random pos.
  x_rand <- sample(1:lattice_size, 1) 
  y_rand <- sample(1:lattice_size, 1) 
  
  lattice[x_rand, y_rand] <- lattice[x_rand, y_rand] + 2
  if(x_rand > 1){
    lattice[x_rand - 1, y_rand] <- lattice[x_rand - 1, y_rand] - 1
  }
  if(y_rand > 1){
    lattice[x_rand, y_rand - 1] <- lattice[x_rand, y_rand - 1] - 1
  }
  
  
  if(t %% 1e5 == 0){
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

## plotting


table(s_values)
hist(s_values, probability = TRUE)
max(s_values)

table(t_values)
hist(t_values)
max(t_values)

table(l_values)
hist(l_values)  
max(l_values)

plot(density(s_values, from = 1, to = 600), log = 'xy')
plot(density(t_values, bw = .5, from = 1, to = 100), log = 'xy')
plot(density(l_values, bw= 0.3), log = 'xy')




## fitting densities
#library(MASS)
#library(fitdistrplus)


### t
t_dataFrame <- as.data.frame(table(t_values))
t_dataFrame$log10_T <- log10(as.numeric(as.character(t_dataFrame$t_values))) # transformation necessary because these values were "factors" before
t_dataFrame$prob <- t_dataFrame$Freq / length(t_values)
t_dataFrame$log10_prob  <- log10(t_dataFrame$prob)

linearMod_t <- lm(log10_prob ~ log10_T, data=t_dataFrame) # , subset=0:4
summary(linearMod_t)

plot(t_dataFrame$log10_T, t_dataFrame$log10_prob) # , xlim = c(0,4)
abline(linearMod_t)


### s
s_dataFrame <- as.data.frame(table(s_values))
s_dataFrame$log10_s <- log10(as.numeric(as.character(s_dataFrame$s_values))) # transformation necessary because these values were "factors" before
s_dataFrame$prob <- s_dataFrame$Freq / length(s_values)
s_dataFrame$log10_prob  <- log10(s_dataFrame$prob)

linearMod_s <- lm(log10_prob ~ log10_s, data=s_dataFrame) # , subset=0:4
summary(linearMod_s)

plot(s_dataFrame$log10_s, s_dataFrame$log10_prob, type = "l") # , xlim = c(0,4)
abline(linearMod_s)

### l
l_dataFrame <- as.data.frame(table(l_values))
l_dataFrame$l_shifted <- as.numeric(as.character(l_dataFrame$l_values)) + 1  # shift l by one in order to be able to get log values
l_dataFrame$log10_l <- log10(l_dataFrame$l_shifted) # transformation necessary because these values were "factors" before
l_dataFrame$prob <- l_dataFrame$Freq / length(l_values)
l_dataFrame$log10_prob  <- log10(l_dataFrame$prob)

linearMod_l <- lm(log10_prob ~ log10_l, data=l_dataFrame) # , subset=0:4
summary(linearMod_l)

plot(l_dataFrame$log10_l, l_dataFrame$log10_prob, type = "l") # , xlim = c(0,4)
abline(linearMod_l)

