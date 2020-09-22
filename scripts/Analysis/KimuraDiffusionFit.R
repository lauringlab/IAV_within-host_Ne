require(tidyverse)
require(bbmle)
require(glue)
source("scripts/lib/kimuraFunctions.R")
source("./scripts/lib/utils.R")


vars<-read_csv("data/Intrahost_initially_present.csv") %>% filter(within_host_time>0,freq1<0.5)
NSvars<-filter(vars,donor_class=="Nonsynonymous")
Svars<-filter(vars,donor_class=="Synonymous")


write_to_summary("NS iSNV n", nrow(NSvars))
write_to_summary("S iSNV n", nrow(Svars))

kimLLFactory<-function(data,gen_time){
  
  kimLL<-function(Ne=30){
    if(Ne<=0){
      return(1e10)
    }
    kimura_LL(data,Ne,neg=T,gen_time=gen_time)  
  }
return(kimLL)
}

#------------ Model fit generation time 6 ----------------------------

kimuraFit<-mle2(kimLLFactory(vars,6),method="L-BFGS-B",lower=c(0))
kimuraFitProf<-profile(kimuraFit)
write_to_summary("Diffusion model 6 Ne:",glue('{kimuraFit@coef} :[{confint(kimuraFitProf)[1]} - {confint(kimuraFitProf)[2]}]'))

kimuraFitNS<-mle2(kimLLFactory(NSvars,6),method="L-BFGS-B",lower=c(0))
kimuraFitNSProf<-profile(kimuraFitNS)
write_to_summary("NS Ne:",glue('{kimuraFitNS@coef} :[{confint(kimuraFitNSProf)[1]} - {confint(kimuraFitNSProf)[2]}]'))

kimuraFitS<-mle2(kimLLFactory(Svars,6),method="L-BFGS-B",lower=c(0))
kimuraFitSProf<-profile(kimuraFitS)
write_to_summary("S Ne:",glue('{kimuraFitS@coef} :[{confint(kimuraFitSProf)[1]} - {confint(kimuraFitSProf)[2]}]'))

#------------ Model fit generation time 12 ----------------------------

kimura12Fit<-mle2(kimLLFactory(vars,12),method="L-BFGS-B",lower=c(0))
kimuraFit12Prof<-profile(kimuraFit12)
write_to_summary("Diffusion model 12 Ne:",glue('{kimuraFit12@coef} :[{confint(kimuraFit12Prof)[1]} - {confint(kimuraFit12Prof)[2]}]'))

kimuraFit12NS<-mle2(kimLLFactory(NSvars,12),method="L-BFGS-B",lower=c(0))
kimuraFit12NSProf<-profile(kimuraFit12NS)
write_to_summary("NS 12 Ne:",glue('{kimuraFit12NS@coef} :[{confint(kimuraFit12NSProf)[1]} - {confint(kimuraFit12NSProf)[2]}]'))

kimuraFit12S<-mle2(kimLLFactory(Svars,12),method="L-BFGS-B",lower=c(0))
kimuraFit12SProf<-profile(kimuraFit12S)
write_to_summary("S 12 Ne:",glue('{kimuraFit12S@coef} :[{confint(kimuraFit12SProf)[1]} - {confint(kimuraFit12SProf)[2]}]'))





