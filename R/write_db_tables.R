#' Updated db
#'
#' @import dplyr
#' @import DBI
#' @import stringr
#' @import DESeq2
#' @importFrom rtracklayer import
#' @importFrom openxlsx read.xlsx
#' @importFrom magrittr "%>%"
#' @importFrom tibble deframe
conn <- nmdtx:::connect_db()
name <- function(base_path = "/Volumes/beegfs/prj/Niels_Gehring/nmd_transcriptome") {
  dbListTables(conn)

  ## METADATA ####
  metadata <- read.csv(file.path(base_path, "phase2/config/metadata_w_files.csv"))
  metadata %<>%
    filter(!is.na(Replicate)) %>%
    mutate(contrasts = ifelse(str_detect(group, 'Luc', negate = TRUE), group, NA)) %>%
    tidyr::fill(contrasts, .direction = 'up')

  id2group <- metadata %>%
    dplyr::select(CCG_Sample_ID, group)

  group_recode <- metadata %>%
    select(Condition, group) %>%
    filter(!str_detect(Condition, "control")) %>%
    distinct() %>%
    tibble::deframe()

  ## junction support from long reads ####
  has_support <- read.csv(
    file.path(base_path, "result_tx_supported_by_nopore_junctions.csv"),
    row.names = 1
  )

  ## DTU ####
  files <- Sys.glob(file.path(base_path, "phase2/results/dtu/dtu*.xlsx"))
  names(files) <- tools::file_path_sans_ext(basename(files))
  dtu <- lapply(files, openxlsx::read.xlsx)
  # dtu <- lapply(dtu, tibble::rownames_to_column)
  dtu <- bind_rows(dtu, .id = "contrasts")
  dtu <- dtu %>% select(!starts_with("countData"))
  dtu <- dtu %>% select(1:8, starts_with("log2fold_"), "genomicData")
  dtu <- dtu %>% mutate(
    log2fold = coalesce(!!!select(., matches("log2fold_|_control$")))
  )
  dtu <- dtu %>% select(!ends_with("control"))
  dtu <- dtu %>% dplyr::rename(
    transcript_id = featureID, gene_name = groupID
  )
  dtu$contrasts <- gsub(
    dtu$contrasts,
    pattern = "dtu_(.*)_vs_(.*)",
    replacement = "\\1"
  )
  dtu$contrasts <- rlang::exec(recode, !!!group_recode, .x = dtu$contrasts)
  stopifnot(all(dtu$contrasts %in% group_recode))

  ## GTF ####
  # gtf <- rtracklayer::import(
  #   file.path(base_path, "phase2/stringtie_merge/merged_each.gtf"))
  gtf <- rtracklayer::import(
    file.path(base_path, "../DFG_seq_Nanopore/GRCh38_90_SIRV_Set3.gtf")
  )

  gtf <- as.data.frame(gtf) %>%
    dplyr::filter(transcript_id %in% unique(dtu$transcript_id))

  dtu %<>%
    left_join(
      gtf %>% select("transcript_name", "transcript_id", "transcript_biotype"),
      by = c("transcript_id"))

  ## Annotation
  anno <- gtf %>%
    filter(type == "transcript") %>%
    dplyr::select(transcript_id, gene_id, transcript_name, gene_name)

  ## DGE ####
  files <- Sys.glob(file.path(base_path, "dge_results/*.xlsx"))
  names(files) <- tools::file_path_sans_ext(basename(files))
  dge <- lapply(
    files,
    openxlsx::read.xlsx
  )
  names(dge) <- names(files)
  dge <- bind_rows(dge, .id = "contrasts")
  dge <- dge %>%
    dplyr::rename(gene_name = gene, gene_id = SYMBOL)
  dge$contrasts <- gsub(
    dge$contrasts,
    pattern = "(.*)_vs_(.*)",
    replacement = "\\1"
  )
  dge$contrasts <- rlang::exec(recode, !!!group_recode, .x = dge$contrasts)
  stopifnot(all(dge$contrasts %in% group_recode))

  ## Gene counts ####
  gene_counts <- readRDS(file.path(base_path, "phase2/data/dge_dds.RDS"))
  gene_counts <- estimateSizeFactors(gene_counts)
  gene_counts <- counts(gene_counts, normalized = TRUE)
  gene_counts %<>%
    as_tibble(rownames = "gene_id") %>%
    left_join(anno %>% select(gene_name, gene_id), by = "gene_id") %>%
    tidyr::pivot_longer(-c(gene_name, gene_id)) %>%
    mutate(name = str_remove(name, "_[12345]")) %>%
    left_join(metadata %>% select(group, contrasts), by=c("name"="group"))


  ## Transcript counts ####
  tx_counts <- readRDS(file.path(base_path, "phase2/data/tx_counts.RDS")) %>%
    DRIMSeq::counts() %>%
    dplyr::rename(transcript_id = feature_id) %>%
    left_join(anno) %>%
    dplyr::filter(transcript_id %in% unique(dtu$transcript_id))  %>%
    select(-c(gene_id)) %>%
    tidyr::pivot_longer(-c(gene_name, transcript_name, transcript_id)) %>%
    collect() %>%
    mutate(name = str_replace_all(name, "[.]", "-") %>% str_remove(., "_[12345]")) %>%
    group_by(gene_name, name) %>%
    mutate(total = sum(value, na.rm = TRUE)) %>%
    filter(total != 0) %>%
    ungroup() %>%
    mutate(usage = value / total) %>%
    left_join(metadata %>% select(group, contrasts), by=c("name"="group"))

  dbWriteTable(conn, "has_support2", has_support, overwrite = TRUE)
  dbWriteTable(conn, "dtu2", dtu, overwrite = TRUE)
  dbWriteTable(conn, "gtf2", gtf, overwrite = TRUE)
  dbWriteTable(conn, "dge2", dge, overwrite = TRUE)
  dbWriteTable(conn, "anno2", anno, overwrite = TRUE)
  dbWriteTable(conn, "metadata", metadata, overwrite = TRUE)
  dbWriteTable(conn, "gene_counts2", gene_counts, overwrite = TRUE)
  dbWriteTable(conn, "tx_counts2", tx_counts, overwrite = TRUE)

  dbDisconnect(conn)
}

