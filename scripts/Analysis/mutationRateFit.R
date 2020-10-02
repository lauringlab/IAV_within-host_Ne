require(tidyverse)
require(bbmle)
require(glue)
source("./scripts/lib/mutation_functions.R")
source("./scripts/lib/utils.R")

# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
meta<-read_csv("./data/all_meta.sequence_success.csv")

accuracy_stringent<-read_csv("./data/accuracy_stringent.csv")
qual<-read_csv("./data/qual.snv.csv", # The quality iSNV
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character(),
                 dup = col_character()
               )) # read in quality variant calls from all 
# B
# sim.df<-read_csv("./results/simulated_fits.csv")

#-------------------------- prep up variants ---------------------------------- 
meta_quality<-meta %>%filter(snv_qualified==TRUE)
# All else being equal only_one will pick the sample with the higher titer
meta_one<-HIVEr::only_one(meta_quality) %>% mutate(DPI = as.numeric(collect-(onset-1),unit="days")) %>%
  filter(!(SPECID %in% c("HS1530","MH8137","MH8390")), DPI>0) # ,"M54062" probably mixed
qual_min <- subset(qual,freq.var<0.5) 
qual_min<-qual_min %>% filter(SPECID %in% meta_one$SPECID) # Only get 1 sample person

meta_one<-meta_one %>% mutate(variants=map(SPECID,~dfToList(filter(qual_min,SPECID==.x))))

write_to_summary("Samples used:",length((meta_one$ENROLLID)))
write_to_summary("Samples without variants",length(meta_one$ENROLLID)-length(unique(qual_min$ENROLLID)))


#-------------------------- Run mutational Analysis --------------------------- 
# # Estimating $N_e$ and $\mu_e$
# The goal here is to allow for the mutations to enter the system at a given 
# rate so that we can go backwards in time when we need to. Here we are relying 
# on the assumption that no polymorphisms are present at the start of infection
# 
# Here the probability function for $0<f<1$ is given by
# 
# 
# $$
#   g(f,t) = \frac{2 \mu N_e}{f} e^{- \frac{2N_ef}{t}}
# $$
#------------------------------------------------------------------------------


base<-seq(1,9,0.5)
power<- -6  #seq(-7,-5,1)
Ne <- seq(10,90,1)

mu.df<-expand.grid(base,power)
mu = mu.df$Var1*10^mu.df$Var2

parameters<-expand.grid(mu,Ne)
names(parameters)<-c("mu","Ne")
muNeData<-dfToList(meta_one)

output<- parameters %>% rowwise()%>%
  mutate(LL = mut_model_LL(muNeData,mu,Ne,gen_time=6,neg=F,acc=accuracy_stringent)) 

LL<-function(mu=-12,Ne=30){
  if(exp(mu)==0){
    return(1e10)
  }
  ll<-mut_model_LL(muNeData,mu=exp(mu),Ne=Ne,neg=T)
  # print(glue("{exp(mu)},{Ne}: {ll}"))
  return(ll)
}

fit<-mle2(LL,method="L-BFGS-B")
prof<-profile(fit)
output<-rbind(output,tibble(mu=c(exp(fit@coef["mu"])),Ne=c(fit@coef["Ne"]),LL=c(-fit@min)))
write_csv(output,"results/mutationRateFit.csv")

write_to_summary("Mu:",glue('{exp(fit@coef["mu"])} : [{exp(confint(prof)["mu",1])} - {exp(confint(prof)["mu",2])}]'))
write_to_summary("Ne_spec:",glue('{fit@coef["Ne"]} :[{confint(prof)["Ne",1]} - {confint(prof)["Ne",2]}]'))
