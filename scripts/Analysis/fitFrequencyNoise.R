# Title     : Frequency noise
# Objective : Fit an approximation to the binomial distribution to describe the 
# noise in the frequency data.
# Created by: jtmccrone
# Created on: 2020-08-26

#Libraries

require(tidyverse)
require(bbmle)
require(glue)
source("./scripts/lib/utils.R")

# Variants used in the analysis
qual_min<-read_csv("data/qual.snv.csv",
                   # These are the iSNV that were used in the study.
                   col_types = list(
                     ENROLLID= col_character(),
                     SPECID = col_character(),
                     LAURING_ID = col_character(),
                     Id = col_character(),
                     dup =col_character()
                   )) %>% filter(freq.var<0.5) # read in quality variant calls from all
# All variants from samples sequenced twice.
freqs<-read_csv("data/duplicate_sequences.csv")

# now we will select only iSNV used in the analysis
# To do this I will make a SPECID_mutation column linking each iSNV with the isolate
# in which it was found.
quality_mut<-paste0(qual_min$mutation, qual_min$SPECID)

(freqs<-freqs %>% mutate(mut_specid = paste0(mutation,SPECID_original), # here we use SPECID_original to match that used in qual above
                         used= mut_specid %in% quality_mut) %>%
    filter(used==T) %>%
    mutate(meanFreq = map2_dbl(freq1, freq2,~ mean(c(.x,.y)))))
dot_plot.discrete<-freqs%>% # filter so only those used are plotted
  ggplot(aes(x=freq1,y=freq2))+
  geom_point(aes(color=as.factor(floor(log10(gc_ul)))))+
  geom_abline(slope=1,intercept = 0,lty=2) +
  xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 2")
dot_plot.discrete


################################################################################
##                                                                            ##
##                       Normal approximation                                 ##
##                    mean  = np and variance (sigma^2)= np(1-p)              ##  
##                However here we are dealing with frequences so              ##
##          mean = p and variance becoes (p-1)/n (variance is devided by n^2) ##
################################################################################

nLL_norm<-function(n=100){
  freqs %>% mutate(nLL = -1*(dnorm(freq1,meanFreq,sqrt(meanFreq*(1-meanFreq)/n),TRUE) 
                             + dnorm(freq2,meanFreq,sqrt(meanFreq*(1-meanFreq)/n),TRUE) ) ) %>%
    select(nLL)%>% sum()
}
mle2(nLL_norm,method="L-BFGS-B")


################################################################################
##                                                                            ##
##                       Beta approximation                                   ##
##      A beta distribution descibes the posterior distribution of p when p   ##
##  is estimated from a binomial distribution with a uniform beta prior on p  ##                                                                   ##
##       α = μν, β = (1 − μ)ν where v is the sample size or sample size +2     ##                                              ##
################################################################################


nLL_beta<-function(v=100){
  freqs %>% mutate(nLL = -1*(dbeta(freq1,meanFreq*v,(1-meanFreq)*v,log= TRUE) + dbeta(freq2,meanFreq*v,(1-meanFreq)*v,log=TRUE))) %>%
    select(nLL)%>% sum()
}

nLL_beta_no_outlier<-function(v){
  freqs %>% filter(freq2>0.1)%>% mutate(nLL = -1*(dbeta(freq1,meanFreq*v,(1-meanFreq)*v,log= TRUE) + dbeta(freq2,meanFreq*v,(1-meanFreq)*v,log=TRUE))) %>%
    select(nLL)%>% sum()
}


betaFit<-mle2(nLL_beta,method="L-BFGS-B")
betaFit

profileBetaFit<-profile(betaFit)

write_to_summary("v:",glue('{betaFit@coef} : [{confint(profileBetaFit)[1]} - {confint(profileBetaFit)[2]}]'))


## Some plots to help intuition
# The variance in the normal fit is slightly lower than the beta but they are 
# comparable


getVar<-function(p,v){
  a<- p*v
  b<-(1-p)*v
  
  (a*b)/((a+b)^2 *(a+b+1))
}
getNormVar<-function(p,n){
  p*(1-p)/n
}
x<-seq(0.001,0.999,by=0.001)

plot(x,getVar(x,503),type="l")

lines(x,getNormVar(x,603),lty=2)



