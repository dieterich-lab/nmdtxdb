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
#' @importFrom tidyr fill
populate_db <- function() {

  base_path = "/Volumes/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/data/"
  conn <- nmdtx:::connect_db()
  dbListTables(conn)


  metadata <- readRDS(file.path(base_path, 'metadata.RDS'))

  metadata %<>%
    filter(!is.na(Replicate)) %>%
    mutate(contrasts = ifelse(str_detect(group, "Luc", negate = TRUE), group, NA)) %>%
    fill(contrasts, .direction = "up")

  copy_to(conn,
          metadata,
          "metadata",
          temporary = FALSE,
          indexes = list(
            "contrasts",
            "CCG_Sample_ID"
          )
  )

  id2group <- metadata %>%
    dplyr::select(CCG_Sample_ID, group)

  group_recode <- metadata %>%
    select(Condition, group) %>%
    filter(!str_detect(Condition, "control")) %>%
    distinct() %>%
    tibble::deframe()

  ## DTU ####
  files <- Sys.glob(file.path(base_path, "phase2/results/dtu/dtu*.xlsx"))
  names(files) <- tools::file_path_sans_ext(basename(files))
  dtu <- lapply(files, openxlsx::read.xlsx)
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
  gtf <- readRDS(file.path(base_path, 'gtf.RDS'))
  gtf <- as.data.frame(gtf)

  copy_to(conn,
          gtf,
          "gtf",
          temporary = FALSE,
          indexes = list(
            "transcript_id",
            "gene_name",
            "ref_gene_id"
          )
  )

  anno <- readRDS(file.path(base_path, 'tx2gene.RDS')) %>%
    select(-keep) %>%
    mutate(gene_name = coalesce(gene_name, gene_id)) %>%
    group_by(gene_id) %>%
    mutate(transcript_name = rank(transcript_id)) %>%
    mutate(transcript_name = str_c(gene_name, transcript_name, sep = '_'))
    # order by expression?

  copy_to(conn,
          anno,
          "anno",
          temporary = FALSE,
          indexes = list(
            "transcript_id",
            "gene_name",
            "gene_id",
            "transcript_name"
          )
  )

  dtu %<>%
    left_join(
      gtf %>% select("transcript_name", "transcript_id", "transcript_biotype"),
      by = c("transcript_id")
    )

  dge <- readRDS(file.path(base_path, 'dge_results.RDS'))
  dge$contrast <- gsub(
    dge$contrast,
    pattern = "(.*)-vs-(.*)",
    replacement = "\\1"
  )
  # dge$contrasts <- rlang::exec(recode, !!!group_recode, .x = dge$contrasts)
  stopifnot(all(dge$contrasts %in% group_recode))
  copy_to(conn,
          anno,
          "dge",
          temporary = FALSE,
          indexes = list(
            "gene_name",
            "gene_id"
        )
  )


  ## Gene counts ####
  base_path <- "/Volumes/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/data"
  txi.genes <- readRDS(file.path(base_path, "gene_counts.RDS"))

  dds <- DESeqDataSetFromTximport(
    txi = txi.genes,
    colData = metadata |> select('group'),
    design = ~group)

  dds <- DESeq(dds)

  gene_counts <- counts(dds) %>%
    as_tibble(rownames = "gene_id") %>%
    tidyr::pivot_longer(-gene_id) %>%
    mutate(name = str_remove(name, "_[12345]")) %>%
    left_join(metadata %>% select(group, contrasts), by = c("name" = "group"))

  copy_to(conn,
          gene_counts,
          "gene_counts",
          temporary = FALSE,
          indexes = list(
            "gene_id"
          )
  )


  ## Transcript counts ####
  tx_counts <- readRDS(file.path(base_path, "drimseq_data.RDS")) %>%
    DRIMSeq::counts() %>%
    dplyr::rename(transcript_id = feature_id) %>%
    left_join(anno %>% select(-ref_gene_id)) %>%
    collect() %>%
    tidyr::pivot_longer(-c(transcript_id, gene_id, gene_name, transcript_name)) %>%
    mutate(name = str_replace_all(name, "[.]", "-") %>% str_remove(., "_[12345]")) %>%
    group_by(gene_id, name) %>%
    mutate(total = sum(value, na.rm = TRUE)) %>%
    filter(total != 0) %>%
    ungroup() %>%
    mutate(usage = value / total) %>%
    left_join(metadata %>% select(group, contrasts), by = c("name" = "group"))
  # see tx_counts to find a few missing contrasts due to misspell groups

  copy_to(conn,
          tx_counts %>% distinct(),
          "tx_counts",
          temporary = FALSE,
          indexes = list(
            "transcript_id",
            "gene_name",
            "gene_id",
            "transcript_name"
          )
  )

  dbDisconnect(conn)
}
