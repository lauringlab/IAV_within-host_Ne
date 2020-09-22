require(tidyverse)
require(cowplot)
require(directlabels)
require(bbmle)
cbPalette<-wesanderson::wes_palette("Zissou") # Set up figures
# --------------------------------- Functions ---------------------------------
write_to_summary<-function(line_pattern,value){
  file = readLines("./results/results.table.tsv")
  line_pattern_regex = paste0("^",line_pattern)
  line = grep(line_pattern_regex,file)
  file[line] = paste0(line_pattern,"\t",value)
  writeLines(file,"./results/results.table.tsv")
}
plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}

# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
# A
mjutationRateAnalysis<-read_csv("./results/")


# --------------------------------- Figure 1A ----------------------------------
#   Likelihood surface for the joint estimation of the mutation rate and Ne
# ------------------------------------------------------------------------------

top<-output[output$LL==max(output$LL),]

tiles.p<-ggplot(output,aes(x=Ne,y=mu,z=LL))+
  scale_y_log10(breaks = seq(1e-6,1e-5,1e-6))+
  stat_contour(breaks = c(-3774,seq(-3807,-4000,by=-30)),aes(color = ..level..),show.legend = F)

  # geom_point(data=top,aes(x=Ne,y=mu))

tiles.p<-direct.label(tiles.p,list("bottom.pieces"))
tiles.p<-tiles.p +  theme(text=element_text(family = 'Helvetica',size=12))
tiles.p
save_plot("./results/Figures/Figure1A.pdf", tiles.p,
          base_aspect_ratio = 1.3)

embed_fonts("./results/Figures/Figure1A.pdf")

write.csv(rename(output,Log.likelihood=LL),"./results/Figures/Figure1A.csv")




# --------------------------------- Figure 1B ----------------------------------
# 
#-------------------------------------------------------------------------------


save_plot("./results/Figures/Figure1B.pdf", fig_1B,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1B.pdf")
write.csv(sim.df.l,"./results/Figures/Figure1B.csv")


