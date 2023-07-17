#' Updated db
#'
#' @import dplyr
#' @import DBI
#' @import stringr
#' @import DESeq2
#' @import tidyr
#' @importFrom magrittr "%>%"
#' @importFrom tibble deframe
#' @importFrom tidyr fill
populate_db <- function() {
  base_path <- "/Volumes/beegfs/prj/Niels_Gehring/nmd_transcriptome/phaseFinal/data/"
  conn <- nmdtx:::connect_db()
  dbListTables(conn)


  metadata <- readRDS(file.path(base_path, "metadata.RDS"))
  metadata <- as.data.frame(metadata)
  metadata$Knockout <- replace_na(metadata$Knockout, "_NoKO")
  metadata$Knockout <- str_replace(metadata$Knockout, "-", "")
  metadata$Knockdown <- str_c(metadata$Knockdown, "KD")
  metadata$group_old <- metadata$group
  metadata$group <- metadata %>%
    dplyr::select(Knockdown, Knockout) %>%
    unite(group, sep = "") %>%
    pull(group)
  metadata$cellline <- word(metadata$Cell_line, 1)
  metadata$group <- str_glue_data(metadata, "{cellline}{Knockout}_{Knockdown}-KD_{clone}", .na = "")
  metadata$group <- gsub(x = metadata$group, "_$", "")
  metadata$group <- as.factor(metadata$group)

  group2contrast <- structure(list(group = c(
    "HEK_NoKO_SMG5KD-KD_Z023", "HEK_NoKO_LucKD-KD_Z023",
    "HEK_NoKO_SMG6+SMG7KD-KD_Z023", "HEK_NoKO_LucKD-KD_Z023", "HEK_SMG7KO_SMG5KD-KD_Z245",
    "HEK_SMG7KO_LucKD-KD_Z245", "HEK_SMG7KO_SMG6KD-KD_Z245", "HEK_SMG7KO_LucKD-KD_Z245",
    "HEK_SMG7KO_SMG5KD-KD_Z319", "HEK_SMG7KO_LucKD-KD_Z319", "HEK_SMG7KO_SMG6KD-KD_Z319",
    "HEK_SMG7KO_LucKD-KD_Z319", "HeLa_NoKO_SMG6+SMG7KD-KD_Z021",
    "HeLa_NoKO_LucKD-KD_Z021", "MCF7_NoKO_SMG6+SMG7KD-KD", "MCF7_NoKO_LucKD-KD",
    "U2OS_NoKO_SMG6+SMG7KD-KD", "U2OS_NoKO_LucKD-KD"
  ), contrasts = c(
    "HEK_NoKO_SMG5KD-KD_Z023-vs-HEK_NoKO_LucKD-KD_Z023",
    "HEK_NoKO_SMG5KD-KD_Z023-vs-HEK_NoKO_LucKD-KD_Z023", "HEK_NoKO_SMG6+SMG7KD-KD_Z023-vs-HEK_NoKO_LucKD-KD_Z023",
    "HEK_NoKO_SMG6+SMG7KD-KD_Z023-vs-HEK_NoKO_LucKD-KD_Z023", "HEK_SMG7KO_SMG5KD-KD_Z245-vs-HEK_SMG7KO_LucKD-KD_Z245",
    "HEK_SMG7KO_SMG5KD-KD_Z245-vs-HEK_SMG7KO_LucKD-KD_Z245", "HEK_SMG7KO_SMG6KD-KD_Z245-vs-HEK_SMG7KO_LucKD-KD_Z245",
    "HEK_SMG7KO_SMG6KD-KD_Z245-vs-HEK_SMG7KO_LucKD-KD_Z245", "HEK_SMG7KO_SMG5KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319",
    "HEK_SMG7KO_SMG5KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319", "HEK_SMG7KO_SMG6KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319",
    "HEK_SMG7KO_SMG6KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319", "HeLa_NoKO_SMG6+SMG7KD-KD_Z021-vs-HeLa_NoKO_LucKD-KD_Z021",
    "HeLa_NoKO_SMG6+SMG7KD-KD_Z021-vs-HeLa_NoKO_LucKD-KD_Z021", "MCF7_NoKO_SMG6+SMG7KD-KD-vs-MCF7_NoKO_LucKD-KD",
    "MCF7_NoKO_SMG6+SMG7KD-KD-vs-MCF7_NoKO_LucKD-KD", "U2OS_NoKO_SMG6+SMG7KD-KD-vs-U2OS_NoKO_LucKD-KD",
    "U2OS_NoKO_SMG6+SMG7KD-KD-vs-U2OS_NoKO_LucKD-KD"
  )), row.names = c(
    NA,
    -18L
  ), class = "data.frame")

  metadata <- metadata %>%
    select(group, cellline, Knockdown, Knockout, group_old) %>%
    mutate(contrasts = group2contrast$contrast[match(group, group2contrast$group)])

  copy_to(conn,
    metadata,
    "metadata",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "contrasts"
    )
  )

  dte <- readRDS(file.path(base_path, "dte_results.RDS")) %>%
    tibble() %>%
    select(1:8, starts_with("log2fold_")) %>%
    mutate(
      log2fold = coalesce(!!!select(., starts_with("log2fold")))
    ) %>%
    select(-starts_with("log2fold_"))
  dte <- dte %>% rename(contrasts = contrast)

  copy_to(conn,
    dte,
    "dte",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "contrasts",
      "gene_id",
      "transcript_id"
    )
  )

  dge <- readRDS(file.path(base_path, "dge_results.RDS")) %>%
    rename(contrasts = contrast)

  copy_to(conn,
    dge,
    "dge",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "contrasts",
      "gene_name",
      "gene_id"
    )
  )

  ## GTF ####
  gtf <- import(file.path(base_path, "gtf_annotated.gtf"))
  gtf <- as.data.frame(gtf)
  gtf$transcript_id <- gsub(x = gtf$transcript_id, "\\_.*", "")
  gtf <- left_join(gtf, anno, by = "transcript_id")

  copy_to(conn,
    gtf,
    "gtf",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "transcript_id",
      "cds_source"
    )
  )
  gff_compare <- read.table(
    file.path(base_path, "../stringtie_merge/fix_comp_ref.merged_each.fix.gtf.tmap"),
    header = 1
  )

  gff_compare %<>%
    select(ref_gene_id, ref_id, qry_gene_id, qry_id) %>%
    rename(transcript_id = ref_id)

  ref_anno <- import("/Volumes/beegfs/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")

  ref_anno <- ref_anno %>%
    as.data.frame() %>%
    filter(type == "transcript") %>%
    select(gene_id, gene_name, transcript_id, transcript_name) %>%
    rename(ref_gene_id = gene_id, ref_gene_name = gene_name, ref_transcript_name = transcript_name)

  gff_compare %<>% left_join(ref_anno)
  gff_compare %<>%
    rename(ref_transcript_id = transcript_id, transcript_id = qry_id, gene_id = qry_gene_id)

  anno <- readRDS(file.path(base_path, "tx2gene.RDS")) %>%
    select(-keep)

  anno$ref_gene_id <- NULL
  anno %<>%
    left_join(gff_compare)

  copy_to(conn,
    anno,
    overwrite = TRUE,
    "anno",
    temporary = FALSE,
    indexes = list(
      "transcript_id",
      "gene_name",
      "gene_id",
      "ref_transcript_id",
      "ref_transcript_name",
      "ref_gene_name",
    )
  )

  dds <- readRDS(file.path(base_path, "gene_counts.RDS")) %>%
    DESeqDataSetFromTximport(
      txi = .,
      colData = metadata |> select("group"),
      design = ~group
    ) %>%
    DESeq(.)

  gene_counts <- counts(dds) %>%
    as_tibble(rownames = "gene_id") %>%
    tidyr::pivot_longer(-gene_id) %>%
    mutate(group = str_remove(name, "_[12345]")) %>%
    left_join(metadata %>% select(group_old, contrasts), by = c("group" = "group_old")) %>%
    distinct()

  copy_to(conn,
    gene_counts,
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "gene_id",
      "contrasts"
    )
  )

  tx_counts <- readRDS(file.path(base_path, "drimseq_data.RDS")) %>%
    DRIMSeq::counts() %>%
    dplyr::rename(transcript_id = feature_id) %>%
    tidyr::pivot_longer(-c(transcript_id, gene_id)) %>%
    mutate(group = str_replace_all(name, "[.]", "-") %>% str_remove(., "_[12345]")) %>%
    group_by(gene_id, group) %>%
    mutate(total = sum(value, na.rm = TRUE)) %>%
    filter(total != 0) %>%
    ungroup() %>%
    mutate(usage = value / total) %>%
    left_join(metadata %>% select(group_old, contrasts), by = c("group" = "group_old")) %>%
    distinct()

  copy_to(conn,
    tx_counts,
    "tx_counts",
    overwrite = TRUE,
    temporary = FALSE,
    indexes = list(
      "transcript_id",
      "gene_id",
      "contrasts"
    )
  )



  dbDisconnect(conn)


}
