# The effective population size and mutation rate of influenza A acutely infected individuals
This repository holds the analysis used in *citation,2020* it relies in part on the [HIVEr R package] (https://github.com/jtmccr1/HIVEr). The code for the BEAST plugin (results/BEAST/plugins/*jar)  is held in a separate repository (https://github.com/jtmccr1/kimuraDiffusionPlugin). Code for the final figures can be found on observablehq https://observablehq.com/collection/@jtmccr1/ne-paper-figures. All other analysis were run in R.  

## Overview
--------
- data : data files used in the analysisas present in commit c0f9740 from Host_Level_IAV_evolution
- scripts: 
  - Analysis: R scripts used to run maximum likelihood models and set up xmls
  - Figures: Scripts used to trial figures
  - lib: R scripts with functions used in Analysis.
- results: Intermediate and final results used to make figures/table
  - BEAST: xml and logs of BEAST analysis
  - results.table.tsv : useful list of maximum likelihood outputs. Numbers may change slightly on reanalysis due to variantion in the optimizer.
  - subsampleiSNVs.RData: RData file with subsampled data set.
  -SubSampleRun_2020-09-28.RData: Output file from running "scripts/Analysis/SubsampleFit.R"
- xmls: holds xml files used in BEAST plugin analysis.  

--------







