## Dev for individual visualizations
suppressPackageStartupMessages({
  library(tidyverse)
  library(reactable)
  library(plotly)
})

conn <- nmdtx:::connect_db()

anno <- function() {
  conn %>% tbl("anno") %>% filter(gene_name == "SRSF2") %>% collect()
}

library(tidyverse)
library(plotly)
log10_or_max <- nmdtx:::log10_or_max


## test components:
library(htmltools)
print(nmdtx:::render_gene_card("ENSG00000124193", conn), browse = T)

profvis::profvis({print(nmdtx::run_app())})

# study reactivity log
library(reactlog)
reactlog_enable()
shiny::reactlogShow()
