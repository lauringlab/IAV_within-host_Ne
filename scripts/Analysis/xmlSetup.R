require(xml2)
require(tidyverse)
require(glue)
meta<-read_csv("./data/all_meta.sequence_success.csv")


intraHostVariants <-read_csv("data/Intrahost_initially_present.csv") %>% filter(freq1<0.5,within_host_time>0)

meta<-meta %>% 
  filter(snv_qualified==TRUE,SPECID %in% c(intraHostVariants$SPECID1,intraHostVariants$SPECID2))%>% 
  mutate(DPI = as.numeric(collect-(onset-1),unit="days"),
         floorLogGcUL = floor(log10(gc_ul)))


doc<-xml_new_root("Beast")
# make Samples
makeSample<-function(data){
  sample<-xml_add_child(doc,"sample")
  xml_attr(sample,"id")<-glue("sample{data$SPECID}")
  xml_attr(sample,"specID")<-data$SPECID
  xml_attr(sample,"DPI")<-data$DPI
  xml_attr(sample,"floorLogGcUL")<-data$floorLogGcUL
  
  variants<-intraHostVariants %>% filter(SPECID1==data$SPECID)
  frequencykey<-"freq1"
  if(nrow(variants)==0){
    variants<-intraHostVariants %>% filter(SPECID2==data$SPECID)
    frequencykey<-"freq2"
    
  }
  variants %>% rowwise() %>%
    group_walk(~(function(x){
      var<-xml_add_child(sample,"isnv")
      xml_attr(var,"id")<-glue('{data$SPECID}_{x$mutation}')
      xml_attr(var,"freq")<-x[frequencykey]
    })(.x))
  
  
}
meta %>% rowwise(SPECID) %>%
  group_walk(~makeSample(.x))

# Make isnv list of initial Frequencies


initialFreqs<-xml_add_child(doc,"isnvList")
xml_attr(initialFreqs,"id")<-"initialFreqs"

intraHostVariants %>% rowwise() %>%
  group_walk(~(function(x){
    var<-xml_add_child(initialFreqs,"isnv")
    xml_attr(var,"idref")<-glue('{x$SPECID1}_{x$mutation}')
  })(.x))


# Make iSNV list of final frequencies
finalFreqs<-xml_add_child(doc,"isnvList")
xml_attr(isnvTraces,"id")<-"finalFreqs"

intraHostVariants %>% rowwise() %>%
  group_walk(~(function(x){
    var<-xml_add_child(initialFreqs,"isnv")
    xml_attr(var,"idref")<-glue('{x$SPECID2}_{x$mutation}')
  })(.x))



# make isnv traces 

isnvTraces<-xml_add_child(doc,"traceList")
xml_attr(isnvTraces,"id")<-"isnvTraces"

intraHostVariants %>% rowwise() %>%
  group_walk(~(function(x){
    trace<-xml_add_child(isnvTraces,"isnvTrace")
    xml_attr(trace,"id")
    obs1<-xml_add_child(trace,"isnv")
    xml_attr(obs1,"idref")<-glue('{x$SPECID1}_{x$mutation}')
    obs2<-xml_add_child(trace,"isnv")
    xml_attr(obs2,"idref")<-glue('{x$SPECID2}_{x$mutation}')

  })(.x))



write_xml(doc,"./test.xml")
