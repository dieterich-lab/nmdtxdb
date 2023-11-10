## Dev for individual visualizations
suppressPackageStartupMessages({
  library(tidyverse)
  library(reactable)
  library(plotly)
})

db <- readRDS("database2.RDS")

anno <- db[['anno']] %>% filter(gene_name == "SRSF2")
contrast <- c("HEK_SMG7KO_SMG6KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319")
tx <- anno$transcript_id


from_start <- readRDS('from_start.RDS')
from_start$thick <- as.character(from_start$thick)
from_start$cdna_thick %>% head()
from_start$is_ptc
from_start[c("blocks",  "cdna_thick.start", "cdna_thick.end", "cdna_thick.width", "itemRgb")] <- NUL




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
