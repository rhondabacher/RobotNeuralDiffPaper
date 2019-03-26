# Scripts for the minutely robot RNA-seq manuscript


This repository contains 

1. Code to reproduce the analyses and figures in "Automated High Frequency Sampling of Stem Cell Differentiation Reveals Early Stage Divergence of Human and Mouse Gene Expression Kinetics". This includes all preprocessing steps and the generation of figures in the main text and supplement. 

2. Trendy output files that can be uploaded to the Trendy shiny application and instructions page.

[Analysis was done using the Trendy R package, available at Bioconductor]: https://bioconductor.org/packages/release/bioc/html/Trendy.html


# Instructions to view data in Trendy Shiny app

### For instructions with screenshots see this page: https://github.com/rhondabacher/RobotSeq/blob/master/TrendyShinyInstructions.pdf


1. Click on the folder "RDATA" above. Then click on the dataset you want to view. Then click "Download".

2. Assuming you have R installed already, type in:

library(BiocManager)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("Trendy")

3. To run the Trendy Shiny App type in:

library(Trendy)
trendyShiny()

4. A webpage will open and will ask you to select the file you want to upload. 



# License

Licensed under GPL-3.
