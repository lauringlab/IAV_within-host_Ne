require(hypergeo)
require(tidyverse)
accuracy_stringent = readr::read_csv("data/accuracy_stringent.csv")
hyp2f1<-function(a,b,c,x){
  o<-hypergeo(a,b,c,x)
  # if(Im(o)==0){
  #   return(Re(o))
  # }
  return(Re(o))
}
ith_term<-function(i,p,t,x,N){ # proofed JT 5/22/17
  q=1-p
  first = p*q*i*(i+1)*((2*i)+1)
  geometric_1= hyp2f1(1-i,i+2,2,p) # verify this the correct function (it's the only hypergeometric I could find with 4 variables) - This is correct
  geometric_2= hyp2f1(1-i,i+2,2,x)
  exponent= i*(i+1)*t/(2*N) # 4N in the book
  out=first*geometric_1*geometric_2*exp(-1*exponent)
  return(out)
}
non_fixed_basic<-function(p,x,t,N,sensitivity=F,gc_ul = NA,acc=NA){
  term = ith_term(i=1,p=p,x=x,t=t,N=N);
  value = term
  i = 2
  while(  i<1000 ){
    term = ith_term(i=i,p=p,x=x,t=t,N=N)
    value= value+term
    
    if((abs(term)<1e-10 & value>0)){
      break
    }
    i=i+1
  }
  if(sensitivity == F){
    return(value)
  }
    # this is the probability of the variant being found where it was given the sensitivity.
    # This is not used in the fitting of the model. All of these variants are found. For each N that we try this term doesn't change. So it is a constant not dependent on N and so doesn't affect the estimate.
    # It is useful though in getting a pdf that sums to 1 in the plots. If we don't have perfect sensitivity for the lost variants then we should treat these the same.
  acc_gc=10**(floor(log10(gc_ul))) # round down to the nearest log10
  if(acc_gc>1e5){
    acc_gc=1e5
    }: # set for the max. We assume we do not gain sensitivity above this cut off. We probably do, but we don't have data on it so this is more conservative
  
  ## Here we assume the accuracy of each range is the same as the smaller range
  if(x<0.05 & x>0.02){
  sense=acc[(acc.gc_ul==acc_gc) & (acc.freq==0.02),'sensitivity']
  sense = sense.iloc[0]
  prob_detect = value*sense
  }else if( x<0.1 & x>0.05){
    sense=acc[(acc.gc_ul==acc_gc) & (acc.freq==0.05),'sensitivity']
  sense = sense.iloc[0]
  prob_detect = value*sense
  }else{
    prob_detect = value # We assume perfect detection above 10%
    
  } 
  
  return(prob_detect)
}

non_fixed<-Vectorize(non_fixed_basic,c("x"))
ith_term_fixed<-function(i,p,t,N){
  # proofed JT 5/22/17
    first = (2*i+1)*p*(1-p)*(-1)**i
  geometric =hyp2f1(1-i,i+2,2,p)
  exponent= i*(i+1)*t/(2*N) # 4N in the book
  out = first*geometric*exp(-1*exponent)
  return(out)
}
below_cut<-function(p,t,N){
    return(integrate(function(x) non_fixed(p,x,t,N),0,0.02)$value) # proofed JT 5/22/17
}
  
just_missed<-function(p,t,N,gc_ul,acc){ # This accounts for the variants that are present but we don't detect them
    acc_gc=10**(floor(log(gc_ul))) # again round down to the nearest log10
  if(acc_gc>1e5){ # set for the max
    acc_gc=1e5
  }
  uncert_term=c()
  f=c(0.02,0.05,0.10)
  for(i in 1:2){ 
    sense=acc$sensitivity[(acc$gc_ul==acc_gc) & (acc$freq==f[i])]
    uncert=1-sense
    # The prob the variant is missed because it is between f[i] and f[i+1] given the sample size
    uncert_term[i]<-integrate(function(x) non_fixed(p,x,t,N),lower=f[i],upper=f[i+1])$value*uncert
  }
  
  return(sum(uncert_term))
}
  
boundaries<-function(p,t,N,final,gc_ul=1e5,sensitivity=F,acc=accuracy_stringent){

  #if final !=0 or final !=1:
  #    raise(ValueError,"Please select 0 or 1 as final frequency")
  if(final==0){
    fixed_freq=1-p   # set for loss. The probabilty the other allele is fixed
  }else if(final ==1){
    fixed_freq = p    # In this case this is the frequency of the allele we want to fix
  }
    
  term = fixed_freq + ith_term_fixed(i=1,p=fixed_freq,t=t,N=N)
  value = term
  i = 2
  while( i<1000 ){
    term = ith_term_fixed(i=i,p=fixed_freq,t=t,N=N)
    value= value+term
    if(abs(term)<1e-10 & value>0){
      break
    }
    i=i+1
  }
  
  if(i==100){
    warning("WOW hit max iterations")
  }

  #print(fixed)
  if(sensitivity == F){
    return(value)
  }
    if(final ==0){
    lost_p = p # this is the frequency of the variant we want to observe loss
    }else if(final==1){
    lost_p = 1-p # in this case we want the loss of the other allele
  }
  below_threshold = below_cut(p=lost_p,t=t,N=N)
  missed= just_missed(p=lost_p,t=t,N=N,gc_ul=gc_ul,acc=acc)
  lost = below_threshold+missed
  return(lost+value)
}
#TODO fix this sensitivity flag
pdf<-function(p,x,t,N,gc_ul=1e5,sensitivity = F, acc=accuracy_stringent){
    if(x <1 & x>0){
    return(non_fixed(p=p,x=x,t=t,N=N,sensitivity=sensitivity,gc_ul=gc_ul,acc=acc))
    }else{
      return(boundaries(p=p,final=x,N=N,t=t,sensitivity=T,gc_ul=gc_ul,acc=acc))
    }
}

kimura_LL<-function(data,Ne,gen_time=6,neg=F,acc=accuracy_stringent){
  dataWithGenerations<-data %>% mutate(generations = within_host_time*(24/gen_time)) %>% dfToList()
  
  LL = map_dbl(dataWithGenerations,log_pdfWrapper,Ne)%>% sum()
  sign= 1
  if(neg){
    sign = -1
  }
  return(sign*LL)
}

log_pdfWrapper<-function(row,Ne){
  return(log(pdf(row$freq1,row$freq2,row$generations,Ne,row$gc_ul2,F)))
}
  
  