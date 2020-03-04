
## plotting


table(simulation_data$s)
hist(simulation_data$s, probability = TRUE)
max(simulation_data$s)

table(simulation_data$t)
hist(simulation_data$t)
max(simulation_data$t)

table(simulation_data$l)
hist(simulation_data$l)  
max(simulation_data$l)

plot(density(simulation_data$s, from = 1, to = 600), log = 'xy')
plot(density(simulation_data$t, bw = .5, from = 1, to = 100), log = 'xy')
plot(density(simulation_data$l, bw= 0.3), log = 'xy')




## fitting densities

### t
t_dataFrame <- as.data.frame(table(simulation_data$t))
t_dataFrame$log10_T <- log10(as.numeric(as.character(t_dataFrame$Var1))) # transformation necessary because these values were "factors" before
t_dataFrame$prob <- t_dataFrame$Freq / sum(t_dataFrame$Freq)
t_dataFrame$log10_prob  <- log10(t_dataFrame$prob)

linearMod_t <- lm(log10_prob ~ log10_T, data=t_dataFrame, subset=(log10_T < 1.9))
summary(linearMod_t)

plot(t_dataFrame$log10_T, t_dataFrame$log10_prob, type = 'l') # , xlim = c(0,4)
abline(linearMod_t)


### s
s_dataFrame <- as.data.frame(table(simulation_data$s))
s_dataFrame$log10_s <- log10(as.numeric(as.character(s_dataFrame$Var1))) # transformation necessary because these values were "factors" before
s_dataFrame$prob <- s_dataFrame$Freq / sum(s_dataFrame$Freq)
s_dataFrame$log10_prob  <- log10(s_dataFrame$prob)

linearMod_s <- lm(log10_prob ~ log10_s, data=s_dataFrame, subset=(log10_s < 2.3)) # , subset=0:4
summary(linearMod_s)

plot(s_dataFrame$log10_s, s_dataFrame$log10_prob, type = "l") # , xlim = c(0,4)
abline(linearMod_s)

### l
l_dataFrame <- as.data.frame(table(simulation_data$l))
l_dataFrame$log10_l <- log10(as.numeric(as.character(l_dataFrame$Var1))) # transformation necessary because these values were "factors" before
l_dataFrame$prob <- l_dataFrame$Freq / sum(l_dataFrame$Freq)
l_dataFrame$log10_prob  <- log10(l_dataFrame$prob)

linearMod_l <- lm(log10_prob ~ log10_l, data=l_dataFrame, subset=(log10_l < 1.3)) # , subset=(SIZE>0.8 & SIZE<7)
summary(linearMod_l)

plot(l_dataFrame$log10_l, l_dataFrame$log10_prob, type = "l") # , xlim = c(0,4)
abline(linearMod_l)



#### conditional expectation values

distinct_t  <- sort(unique(simulation_data$t))
distinct_s  <- sort(unique(simulation_data$s))
distinct_l  <- sort(unique(simulation_data$l))


# conditional on s:
cond_on_s <- data.frame(s=numeric(0),E_t_on_s=numeric(0),E_l_on_s=numeric(0))

s_t_table <- table(simulation_data$s, simulation_data$t) # t sind so spalten und s reihen in table
s_l_table <- table(simulation_data$s, simulation_data$l)


for (s_ in as.character(distinct_s)) {
  
  E_t_on_s_ <- sum(distinct_t * s_t_table[s_,]/sum(s_t_table[s_,]))
  E_l_on_s_ <- sum(distinct_l * s_l_table[s_,]/sum(s_l_table[s_,]))
  
  cond_on_s <- rbind(cond_on_s, 
                     data.frame(s=as.numeric(s_), E_t_on_s=E_t_on_s_ ,E_l_on_s=E_l_on_s_))
}


# conditional on t:
cond_on_t <- data.frame(t=numeric(0),E_s_on_t=numeric(0),E_l_on_t=numeric(0))

t_s_table <- table(simulation_data$t, simulation_data$s) 
t_l_table <- table(simulation_data$t, simulation_data$l)


for (t_ in as.character(distinct_t)) {
  
  E_s_on_t_ <- sum(distinct_s * t_s_table[t_,]/sum(t_s_table[t_,]))
  E_l_on_t_ <- sum(distinct_l * t_l_table[t_,]/sum(t_l_table[t_,]))
  
  cond_on_t <- rbind(cond_on_t, 
                     data.frame(t=as.numeric(t_), E_s_on_t=E_s_on_t_ ,E_l_on_t=E_l_on_t_))
}

# conditional on l:
cond_on_l <- data.frame(l=numeric(0),E_s_on_l=numeric(0),E_t_on_l=numeric(0))

l_s_table <- table(simulation_data$l, simulation_data$s) 
l_t_table <- table(simulation_data$l, simulation_data$t)


for (l_ in as.character(distinct_l)) {
  
  E_s_on_l_ <- sum(distinct_s * l_s_table[l_,]/sum(l_s_table[l_,]))
  E_t_on_l_ <- sum(distinct_t * l_t_table[l_,]/sum(l_t_table[l_,]))
  
  cond_on_l <- rbind(cond_on_l, 
                     data.frame(l=as.numeric(l_), E_s_on_l=E_s_on_l_ ,E_t_on_l=E_t_on_l_))
}

###############################################################
######## plotting and fitting of conditional expectation values
###############################################################
##### cond. s:

cond_on_s$log10_s <- log10(cond_on_s$s)
cond_on_s$log10_E_t_on_s <- log10(cond_on_s$E_t_on_s)
cond_on_s$log10_E_l_on_s <- log10(cond_on_s$E_l_on_s)

# E(t | S=s)
linearMod_E_s_1 <- lm(log10_E_t_on_s ~ log10_s, data=cond_on_s, subset=(log10_s < 3.0))
summary(linearMod_E_s_1)
plot(cond_on_s$log10_s, cond_on_s$log10_E_t_on_s, type = "l", col = "green4") 
abline(linearMod_E_s_1)

# E(l | S=s)
linearMod_E_s_2 <- lm(log10_E_l_on_s ~ log10_s, data=cond_on_s
                      , subset=(log10_s > 0.3 & log10_s < 2.8))
summary(linearMod_E_s_2)
plot(cond_on_s$log10_s, cond_on_s$log10_E_l_on_s, col = "green4", type = "l") 
abline(linearMod_E_s_2)

##############
##### cond. t:
cond_on_t$log10_t <- log10(cond_on_t$t)
cond_on_t$log10_E_s_on_t <- log10(cond_on_t$E_s_on_t)
cond_on_t$log10_E_l_on_t <- log10(cond_on_t$E_l_on_t)

# E(s | T=t)
linearMod_E_t_1 <- lm(log10_E_s_on_t ~ log10_t, data=cond_on_t, subset=(log10_t < 2.5))
summary(linearMod_E_t_1)
plot(cond_on_t$log10_t, cond_on_t$log10_E_s_on_t, type = "l", col = "green4")
abline(linearMod_E_t_1)

# E(l | T=t)
linearMod_E_t_2 <- lm(log10_E_l_on_t ~ log10_t, data=cond_on_t
                      , subset=(log10_t > 0.3 & log10_t < 1.8))
summary(linearMod_E_t_2)
plot(cond_on_t$log10_t, cond_on_t$log10_E_l_on_t, col = "green4", type = "l") 
abline(linearMod_E_t_2)

###############
##### cond. l:
cond_on_l$log10_l <- log10(cond_on_l$l)
cond_on_l$log10_E_s_on_l <- log10(cond_on_l$E_s_on_l)
cond_on_l$log10_E_t_on_l <- log10(cond_on_l$E_t_on_l)

# E(s | L=l)
linearMod_E_l_1 <- lm(log10_E_s_on_l ~ log10_l, data=cond_on_l, subset=(log10_l < 1.3))
summary(linearMod_E_l_1)
plot(cond_on_l$log10_l, cond_on_l$log10_E_s_on_l, col = "green4", type = "l") # , xlim = c(0,4)
abline(linearMod_E_l_1)

# E(t | L=l)
linearMod_E_l_2 <- lm(log10_E_t_on_l ~ log10_l, data=cond_on_l
                      , subset=( log10_l < 1.4))
summary(linearMod_E_l_2)
plot(cond_on_l$log10_l, cond_on_l$log10_E_t_on_l, col = "green4")#, type = "l") # , xlim = c(0,4)
abline(linearMod_E_l_2)

