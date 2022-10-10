library(DBI)
library(tidyverse)

base_path <- "/Volumes/beegfs/prj/Niels_Gehring/nmd_transcriptome"

conn <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "nmd_transcriptome",
  host =  Sys.getenv("NMD_PGHOST"),
  port =  Sys.getenv("NMD_PGPORT"),
  password = Sys.getenv("NMD_PGPASSWORD"),
  user = Sys.getenv("NMD_PGUSER"),
)

dbListTables(conn)


## junction support from long reads ####
has_support <- read.csv(
  file.path(base_path, "result_tx_supported_by_nopore_junctions.csv"),
  row.names = 1)

## DTU ####
files <-  Sys.glob(file.path(base_path, "phase2/results/dtu*.xlsx"))
names(files) <- tools::file_path_sans_ext(basename(files))
dtu <- lapply(files, openxlsx::read.xlsx)
#dtu <- lapply(dtu, tibble::rownames_to_column)
dtu <- bind_rows(dtu, .id='contrasts')
dtu <- dtu %>% select(!starts_with("countData"))
dtu <- dtu %>% select(1:8, starts_with('log2fold_'), 'genomicData')
dtu <- dtu %>% mutate(
  log2fold = coalesce(!!! select(., matches("log2fold_|_control$"))))
dtu <- dtu %>% select(!ends_with('control'))
dtu <- dtu %>% dplyr::rename(transcript_id = featureID, gene_name=groupID)

## GTF ####
# gtf <- rtracklayer::import(
#   file.path(base_path, "phase2/stringtie_merge/merged_each.gtf"))

gtf <- rtracklayer::import(
  file.path(base_path, "../DFG_seq_Nanopore/GRCh38_90_SIRV_Set3.gtf"))

gtf <- as.data.frame(gtf) %>%
  dplyr::filter(transcript_id %in% unique(dtu$transcript_id))

## Annotation
anno <- gtf %>%
  filter(type == 'transcript') %>%
  dplyr::select(transcript_id, gene_id, transcript_name, gene_name)

## DGE ####
files <-  Sys.glob(file.path(base_path, "dge_results/*.xlsx"))
names(files) <- tools::file_path_sans_ext(basename(files))
dge <- lapply(
  files,
  openxlsx::read.xlsx)
names(dge) <- names(files)
dge <- bind_rows(dge, .id='contrasts')
dge <- dge %>%
  rename(gene_name=gene, gene_id=SYMBOL)

## METADATA ####
metadata <- read.csv(file.path(base_path, 'phase2/config/metadata_w_files.csv'))
metadata %<>%
  filter(!is.na(Replicate)) %>%
  mutate(
    condition = gsub(x=Condition, " ", ""),
    cellline = word(Cell_line, 1),
    group = str_glue("{condition}_{Replicate}"))

id2group <- metadata %>%
  dplyr::select(CCG_Sample_ID, group)

## Gene counts ####
gene_counts <- readRDS(file.path(base_path, 'phase2/data/dge_dds.RDS'))
gene_counts <- estimateSizeFactors(gene_counts)
gene_counts <- counts(gene_counts, normalized=TRUE)
colnames(gene_counts) <- metadata$group
gene_counts %<>%
  as_tibble() %>%
  mutate(gene_name = rownames(.))  %>%
  left_join(anno[c("gene_name", "gene_id")], by="gene_name")

## Transcript counts ####
tx_counts <- readRDS(file.path(base_path, 'phase2/data/tx_counts.RDS')) %>%
  counts() %>%
  dplyr::rename(transcript_id=feature_id, gene_name=gene_id) %>%
  left_join(anno) %>%
  dplyr::filter(transcript_id %in% unique(dtu$transcript_id))

dbWriteTable(conn, "has_support2", has_support, overwrite=TRUE)
dbWriteTable(conn, "dtu2", dtu, overwrite=TRUE)
dbWriteTable(conn, "gtf2", gtf, overwrite=TRUE)
dbWriteTable(conn, "dge2", dge, overwrite=TRUE)
dbWriteTable(conn, "anno2", anno, overwrite=TRUE)
dbWriteTable(conn, "gene_counts2", gene_counts, overwrite=TRUE)
dbWriteTable(conn, "tx_counts2", tx_counts, overwrite=TRUE)

dbListTables(conn)
dbDisconnect(conn)
