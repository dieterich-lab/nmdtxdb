options(golem.app.prod = FALSE, warn = 1)
options(shiny.port = httpuv::randomPort())
golem::detach_all_attached()
golem::document_and_reload()
run_app()

## Dev for individual vizualizations
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
