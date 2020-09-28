require(tidyverse)
require(bbmle)
require(glue)
require(furrr)
require(progressr)
source("scripts/lib/kimuraFunctions.R")
source("./scripts/lib/utils.R")


vars<-read_csv("data/Intrahost_initially_present.csv") %>% filter(within_host_time>0,freq1<0.5)
vars %>%
  group_by(ENROLLID) %>%   # prep for work by Species
  nest() %>%              # --> one row per Species
  ungroup() %>% mutate(samp = map(data,slice_sample,size=1)) %>%
  select(samp) %>% unnest(samp)

# ------Seting up Analysis ----------------
# make list of 1000 dataframes each with 1 var/person
runs<-list()
for(i in 1:100){
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

plan(multisession, workers = 4)
# ------Running Model ----------------

with_progress({
  p <- progressor(steps = length(runs))
  results<-future_map(runs,~{
    p()
    mle2(kimLLFactory(.x,6),method="L-BFGS-B",lower=c(0))
    },.progress = T)
})

Ne<-map_dbl(results,~coef(.x)[1])

write_to_summary("Subset median 6 (100) Ne",median(Ne))
write_to_summary("Subset IQR 6 (100) Ne",IQR(Ne))

