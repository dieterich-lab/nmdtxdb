## Dev for individual visualizations
suppressPackageStartupMessages({
  library(tidyverse)
  library(reactable)
  library(plotly)
})

# db <- readRDS("database2.RDS")
devtools::load_all()
anno <- db[['anno']] %>% filter(gene_name == "BAG3")
contrast <- c("HEK_SMG7KO_SMG6KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319")
tx <- anno$transcript_id
cds <- c('canonical', 'ensembl')

library(tidyverse)
library(plotly)
log10_or_max <- nmdtx:::log10_or_max


## test components:
library(htmltools)
print(nmdtx:::render_gene_card("ENSG00000124193"), browse = T)

profvis::profvis({print(nmdtx::run_app())})

# study reactivity log
library(reactlog)
reactlog_enable()
shiny::reactlogShow()
