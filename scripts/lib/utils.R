write_to_summary<-function(line_pattern,value){
  file = readLines("./results/results.table.tsv")
  line_pattern_regex = paste0("^",line_pattern)
  line = grep(line_pattern_regex,file)
  file[line] = paste0(line_pattern,"\t",value)
  writeLines(file,"./results/results.table.tsv")
}
dfToList<-function(df){
  return(split(df,seq(nrow(df))))
}
