## length function
prob_la <- function(Lmax = 250,ages = 2:30,binsize = 5,dat = dat) {
  #set lenbins and ages vectors
  lenbins <- seq(1, Lmax, by = binsize)
  #empty array for probabilities
  prob_len_at_age <- array(0,c(2,length(lenbins),length(ages)))
  
  for (a in 1:length(ages)) {
    for (s in 1:2) {
      #index dat
      tmp <- dat[which(dat$AgeClass == a + 1 & dat$sex == s),]
      
      # calculate probabilites for each length bin
      for (l in 1:length(lenbins)) {
        # first length bin
        if(l == 1){
          prob_len_at_age[s,l,a] <- pnorm(lenbins[l],mean = tmp$mean,
                                          sd = tmp$sd)
          # maximum length bin
        } else if(l == length(lenbins)){
          prob_len_at_age[s,l,a] <- 1 - pnorm(lenbins[l],mean = tmp$mean,
                                              sd = tmp$sd)
        }else{
          # all other length bins
          prob_len_at_age[s,l,a] <- pnorm(lenbins[l],mean = tmp$mean,
                                          sd = tmp$sd) - pnorm(lenbins[l - 1],mean = tmp$mean,
                                                               sd = tmp$sd)
        }
      }
    }
  }
  #return prob array
  prob_len_at_age
}
