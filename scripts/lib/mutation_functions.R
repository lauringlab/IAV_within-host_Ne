# TODO rewrite so not so much messaging
#' Mutation rate and Ne from spectrum
#'
#' This function determines the Ne and mutation rate assuming all infections begin as a clone
#' and then mutations are built upon this base. The probability of a site going from 0 to
#' x given mu, Ne, and t generations given in Rouzine 2001 JVI as
#'
#' \deqn{g(f,t) = \frac{2 \mu N_e}{f} e^{- \frac{2N_ef}{t}}}
#'
#' This model makes some assumptions.
#' 1) Every infection began with every potential loci with a mutation at frequency 0 - we only model minor alleles
#' 2) The population develops under neutral processes
#' 3) Because we integrate numerically we assume anything present at less than 0.001 (0.1%)
#' is essentially at 0.
#' 4) There are 13133 sites that can contain mutations - the number of nucleotides in our concatenated 2014-2015 sample
#'
#' @param data The data frame containing iSNV calls. It must have the following columns:
#' freq.var,DPI-days post infection, and only contain 1 allele/ site. Sites not present will be assumed
#' to not have a variant present.
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param gen_time The generation time in hours. Default to 6
#' @param neg Boolean return the negative log likelihood not. Default =FALSE
#' @param acc An accuracy data frame
#' @return The log likelihood for the data given the parameters and the model
#'
#' @export
dplyr.summarise.inform<-F
accuracy_stringent = readr::read_csv("data/accuracy_stringent.csv")

mut_model_LL<-function(data,mu,Ne,gen_time=6,neg=FALSE,acc=accuracy_stringent){
  
  LL=map_dbl(data,calculateSampleLL,mu,Ne,gen_time) %>% sum()
  if(neg==T){
    return(-(LL))
  }else{
    return(LL) # return total likelihood for found and missed samples.
  }
}
calculateSampleLL<-function(data,mu,Ne,gen_time){

  #First we get the contribution of each iSNV that was found
  dpi = data$DPI
  gc_ul = data$gc_ul

  found_LL<-map_dbl(data$variants[[1]],calculateVariantLL,mu,Ne,dpi*(24/gen_time),gc_ul) %>% sum()
 
  sitesWithoutVariants = 13133-length(data$variants[[1]]);
  sample_LL = found_LL + sitesWithoutVariants * log(not_detected(mu, Ne, t = dpi*(24/gen_time),gc_ul = gc_ul))
  return(sample_LL)
  }

calculateVariantLL<-function(var,mu,Ne,t,gc_ul){
  if(nrow(var)==0){
    return(0)
  }
  return(log(g_ft(mu=mu,Ne=Ne,t=t,
       x=var$freq.var,gc_ul = gc_ul,acc=accuracy_stringent)))
}

#' Likelihood of frequency starting at 0
#'
#' This function determines the Ne and mutation rate assuming all infections begin as a clone
#' and then mutations are built upon this base. The probability of a site going from 0 to
#' x given mu, Ne, and t generations given in Rouzine 2001 JVI as
#'
#' \deqn{g(f,t) = \frac{2 \mu N_e}{f} e^{- \frac{2N_ef}{t}}}
#'
#' The output is adjusted to account for the sensitivity. This isn't important in fitting,
#' but it is nice in testing to know the prob of being found+ prob of not being found =1
#' @param x The measured frequency - must be between 0 and 1
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param t Generations
#' @param gc_ul titer default to NULL use if interested in adjusting for sensitivity
#' @param acc sensitivity data frame default is NULL use if interested in adjusting for sensitivity
#' @export



g_ft<-function(x,mu,Ne,t,gc_ul=NULL,acc=NULL){ # x is the frequency
  prob <- (2*mu*Ne/x) * exp(-1*2*Ne*x/t)
  if(!is.null(acc)){
    acc_gc=10^(floor(log10(gc_ul)))
    if(acc_gc>1e5){
      acc_gc <- 1e5
    }
    if(x>=0.02 & x<0.05){
      sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.02)]
    }else if(x>=0.05 & x<0.1){
      sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.05)]
    }
    else {
      sense=1
    }
    return( prob*sense )
  }
  
  else{
    return(prob)
  }
}



#' Mutational model- non polmorphic
#'
#' What's the probability the site still is at frequency 0 i.e. below the level of detection?
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param t Generations

p0 <- function(mu,Ne,t,limit=0.02){ # 1 - polymorphic - Polymorphic being defined as between 0.001 and 1
  poly <- integrate(g_ft,lower=limit,upper=1,mu=mu,Ne=Ne,t=t)
  1-poly$value
}

#' Probability of not detecting a polymorphism
#' This is the probability the site is less than 0.001 (assumed 0 in this model) or
#' not found by our sequencing method.
#'
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param t The number of generations
#' @param gc_ul The sample titer
#' @param acc  dataframe holding the sensitivity
#' @return The probability for the data given the parameters and the model

not_detected<-function(mu,Ne,t,gc_ul,acc=accuracy_stringent){
  # p0
  naught = p0(mu,Ne,t) # probability it is at 0 or below the limit of detection. 
  # # below cut
  # below = integrate(g_ft,lower=0.001,upper=0.02,
  #                   mu=mu,Ne=Ne,t=t) # below detection
  # missed
  acc_gc=10^(floor(log10(gc_ul)))
  if(acc_gc>1e5){
    acc_gc <- 1e5
  } # This will round the gc down to the nearest log as we have always done to be conservative
  uncert_term=c()
  f=c(0.02,0.05,0.10)
  for(i in 1:(length(f)-1)){
    uncert=1-acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==f[i])]
    # The prob the variant is missed because it is between f[i] and f[i+1] given the sample size
    integrand = integrate(g_ft,lower=f[i],upper=f[i+1],
                          mu=mu,Ne=Ne,t=t)
    uncert_term[i]=integrand$value*uncert # the probability it is present * the probability of not seeing it
  }
  missed = sum(uncert_term)
  return(naught+missed)
}





