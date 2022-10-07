library(DBI)

base_path <- "/Volumes/beegfs/prj/Niels_Gehring/nmd_transcriptome"

conn <- dbConnect(
  RPostgres::Postgres(),
  dbname = "nmd_transcriptome",
  host = "***REMOVED***",
  port = ***REMOVED***,
  password = "***REMOVED***",
  user = "***REMOVED***"
)

dbListTables(conn)

## DTU ####

dtu <- read_csv("~/repos/NMD-Transcriptome/results/phase1/DEXSeq_dtu.csv.gz")
dtu <- dtu %>% group_split(contrasts)
all(dtu[[1]]$transcript_id == dtu[[2]]$transcript_id)
cols <- c("SMG6kd_SMG7ko", "log2fold_SMG6kd_SMG7ko_control", "countData.7", "countData.8", "countData.9")
dtu <- cbind(dtu[[1]] %>% dplyr::select(-all_of(cols)), dtu[[2]][cols])

## GTF ####
gtf <- rtracklayer::import("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf")
gtf <- as.data.frame(gtf) %>%
  dplyr::filter(transcript_id %in% unique(dtu$transcript_id))


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
metadata <- readRDS(here('data/metadata.RDS'))
metadata[1, 1] <- "A006200072_119952"
metadata <- filter(
  metadata, VB_Internal_Sample_ID %in%
    c("33F-1", "33F-2", "33F-3",
      "33F-10", "33F-11", "33F-12",
      "33F-13", "33F-14", "33F-15")) %>%
  mutate(group = str_glue("{group}_{Replicate}"))

id2group <- metadata %>%
  dplyr::select(CCG_Sample_ID, group) %>%
  deframe()

## Gene counts ####
gene_counts <- readRDS(here('phase1/data/salmon_counts_genes.RDS'))$counts %>%
  as_tibble(rownames = "gene_id") %>%
  dplyr::select(c("gene_id", metadata$CCG_Sample_ID)) %>%
  dplyr::filter(gene_id %in% unique(dtu$gene_id)) %>%
  rename(id2group)

## Transcript counts ####
tx_counts <- readRDS(here('phase1/data/salmon_counts.RDS'))$counts %>%
  as_tibble(rownames = "transcript_id") %>%
  left_join(., anno[c('gene_id', 'transcript_id')]) %>%
  dplyr::select(c("gene_id", "transcript_id", metadata$CCG_Sample_ID)) %>%
  dplyr::filter(transcript_id %in% unique(dtu$transcript_id)) %>%
  rename(id2group)

## Annotation
anno <- gtf %>%
  filter(type == 'transcript') %>%
  dplyr::select(transcript_id, gene_id, transcript_name, gene_name)


dbWriteTable(conn, "dtu", dtu, overwrite=TRUE)
dbWriteTable(conn, "gtf", gtf, overwrite=TRUE)
dbWriteTable(conn, "dge", dge, overwrite=TRUE)
dbWriteTable(conn, "anno", anno, overwrite=TRUE)
dbWriteTable(conn, "gene_counts", gene_counts, overwrite=TRUE)
dbWriteTable(conn, "tx_counts", tx_counts, overwrite=TRUE)

dbListTables(conn)
dbDisconnect(conn)
