require(tidyverse)
require(cowplot)
require(directlabels)
require(bbmle)
require(r2d3)
# cbPalette<-wesanderson::wes_palette("Zissou") # Set up figures
# --------------------------------- Functions ---------------------------------
# write_to_summary<-function(line_pattern,value){
#   file = readLines("./results/results.table.tsv")
#   line_pattern_regex = paste0("^",line_pattern)
#   line = grep(line_pattern_regex,file)
#   file[line] = paste0(line_pattern,"\t",value)
#   writeLines(file,"./results/results.table.tsv")
# }
# plot.median <- function(x) {
#   m <- median(x)
#   c(y = m, ymin = m, ymax = m)
# }

# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
# A
mutationRateAnalysis<-read_csv("./results/mutationRateFit.csv")


# --------------------------------- Figure 1A ----------------------------------
#   Likelihood surface for the joint estimation of the mutation rate and Ne
# ------------------------------------------------------------------------------

top<-mutationRateAnalysis[mutationRateAnalysis$LL==max(mutationRateAnalysis$LL),]

tiles.p<-ggplot(mutationRateAnalysis,aes(x=Ne,y=mu,z=LL))+
  scale_y_log10(breaks = seq(1e-6,1e-5,1e-6))+
  stat_contour(breaks = c(-3315,seq(-3320,-4000,by=-50)),aes(color = ..level..),show.legend = F) +
  geom_point(data=top,aes(x=Ne,y=mu)) + theme_cowplot()

# tiles.p<-direct.label(tiles.p,method=list("bottom.pieces"))
tiles.p<-tiles.p +  theme(text=element_text(family = 'Helvetica',size=12))
tiles.p
save_plot("./results/Figures/Figure1A.pdf", tiles.p,
          base_aspect_ratio = 1.3)

embed_fonts("./results/Figures/Figure1A.pdf")

# write.csv(rename(output,Log.likelihood=LL),"./results/Figures/Figure1A.csv")




# --------------------------------- Figure 1B ----------------------------------
# 
#-------------------------------------------------------------------------------

hpd<-c(33,72)


kderaw<-read_tsv("results/initial.extinct.Ne.kde.tsv") %>% rename(Ne="combined.initial.extinct.1M.log",posterior = "NeParameter") 

kde<-kderaw %>%
  # rbind(tibble(Ne=seq(0,min(kderaw$Ne)),posterior=0)) %>% 
  # rbind(tibble(Ne=seq(max(kderaw$Ne),100),posterior=0)) %>% # fill in missing values
  mutate(prior = dgamma(Ne,0.036,scale=1000)) %>% 
  pivot_longer(!Ne, names_to = "type", values_to = "density") %>% 
  mutate(
    fill = case_when(
      # type== "prior" ~ "prior",
      (type=="posterior" & Ne>hpd[1] &Ne<hpd[2]) ~"posterior",
    )
    )
  

r2d3(data=kde, script = "scripts/Figures/Figure1x.js")





fig_1B<-ggplot(data=kde)+geom_line(aes(x=Ne,y=density,color=type),size=1) + 
  geom_ribbon(aes(ymax=density,ymin=0,x=Ne,fill=fill),alpha=0.3) + 
  scale_color_brewer(type="qual","") + scale_fill_brewer(type="qual","") + theme_cowplot()
fig_1B
save_plot("./results/Figures/Figure1B.pdf", fig_1B,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1B.pdf")


