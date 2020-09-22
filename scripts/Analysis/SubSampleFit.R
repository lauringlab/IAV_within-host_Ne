require(tidyverse)
require(bbmle)
require(glue)
require(furrr)
source("scripts/lib/kimuraFunctions.R")
source("./scripts/lib/utils.R")


vars<-read_csv("data/Intrahost_initially_present.csv") %>% filter(within_host_time>0,freq1<0.5)
vars %>%
  group_by(ENROLLID) %>%   # prep for work by Species
  nest() %>%              # --> one row per Species
  ungroup() %>% mutate(samp = map(data,slice_sample,size=1)) %>%
  select(samp) %>% unnest(samp)

# make list of 1000 dataframes each with 1 var/person
runs<-list()
for(i in 1:1000){
  runs[[i]]<-vars %>%
    group_by(ENROLLID) %>%   # prep for work by Species
    nest() %>%              # --> one row per Species
    ungroup() %>% mutate(samp = map(data,slice_sample,size=1)) %>%
    select(samp) %>% unnest(samp)
}


kimLLFactory<-function(data,gen_time){
  
  kimLL<-function(Ne=30){
    if(Ne<=0){
      return(1e10)
    }
    kimura_LL(data,Ne,neg=T,gen_time=gen_time)  
  }
  return(kimLL)
}

plan(multisession, workers = 10)

results<-future_map(runs,~mle2(kimLLFactory(.x,6),method="L-BFGS-B",lower=c(0)))

Ne<-map_dbl(results,~coef(.x)[1])

write_to_summary("Subset median 6 Ne:",median(Ne))
write_to_summary("Subset IQR 6 Ne:",IQR(Ne))

