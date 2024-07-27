## length function
quant_la <- function(A = 30,nquant = 4,dat = dat) {
  #set lenbins and ages vectors
  quants <- seq(0.01, 0.99, length.out = nquant)
  ages <- seq(2,A)
  #empty array for probabilities
  quant_len_at_age <- array(0,c(length(ages),nquant,2))
  
  for (s1 in 1:2) {
    for (a1 in 1:length(ages)) {
      #index dat
      tmp <- dat[which(dat$AgeClass == a1 + 1 & dat$sex == s1),]
      
      quant_len_at_age[a1,,s1] <- qnorm(quants,mean = tmp$mean,sd = tmp$sd)
    }
  }
  
  quant_len_at_age
}
