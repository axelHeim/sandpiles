
RANDU <- function(old_X_i , range){
  X_i <- (old_X_i * 65539) %% (2**31)
  random_number <- ceiling(X_i/(2**31 - 1) * range)
  
  return (c(random_number, X_i))
}

seed <- 1


arr <- c()
for (i in 1:1e6) {
  tmp <- RANDU(seed, 100)
  seed <- tmp[2]
  arr <- c(arr, tmp[1])
}
hist(arr, breaks = 10)
min(arr)
max(arr)

x <- sample(1:100, 1e6, replace = T)
hist(x, breaks = 10)
