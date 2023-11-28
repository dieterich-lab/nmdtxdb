library(ggplot2)
library(patchwork)
library(extrafont)
library(scales)

theme_Publication <- function(base_size = 14, base_family = "Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
  + theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2)
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

plot_figures <- function(db = readRDS("database.RDS")) {
  transcripts <- db[["gtf"]] %>%
    filter(type == "transcript", cds_source != "transdecoder") %>%
    mutate(PTC = color == "#FF0000") %>%
    select(transcript_id, PTC) %>%
    group_by(transcript_id) %>%
    summarise(PTC = any(PTC))

  df <- db[["anno"]] %>%
    mutate(
      class_code = case_when(
        class_code == "=" ~ "Complete",
        class_code == "n" ~ "Retained introns",
        class_code %in% c("c", "j", "k") ~ "Splicing variants",
        TRUE ~ "Other"
      )
    )

  x <- left_join(df, transcripts, by = "transcript_id") %>%
    count(class_code, PTC)


  overview <- x %>%
    mutate(class_code = forcats::fct_reorder(class_code, desc(n))) %>%
    ggplot(aes(fill = PTC, y = n, x = class_code, label = n)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_text(size = 5, position = position_dodge(0.9), color = "black", vjust = -0.2) +
    scale_fill_manual(
      values = c("black", "firebrick", "gray50"),
      name = "is_PTC"
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      n.breaks = 5L,
      expand = c(0, 0, 0.1, 0)
    ) +
    theme_linedraw() +
    theme_Publication() +
    labs(
      title = "NMDtxDB",
      subtitle = "Dataset overview (excluding TransDecoder)",
      x = "Reference match",
      y = "Numb. of transcripts"
    )

  ggsave(
    filename = "~/Nextcloud2/nmd_transcriptome_figures/overview.pdf",
    plot = overview,
    dpi = "retina"
  )

  transcripts <- db[["gtf"]] %>%
    filter(type == "transcript") %>%
    mutate(PTC = color == "#FF0000") %>%
    select(transcript_id, PTC, cds_source)

  df <- db[["anno"]] %>%
    mutate(
      class_code = case_when(
        class_code == "=" ~ "Complete",
        class_code == "n" ~ "Retained introns",
        class_code %in% c("c", "j", "k") ~ "Splicing variants",
        TRUE ~ "Other"
      )
    )

  x <- left_join(df, transcripts, by = "transcript_id", multiple = "all") %>%
    count(class_code, PTC, cds_source) %>%
    distinct()

  overview2 <- x %>%
    mutate(class_code = forcats::fct_reorder(class_code, desc(n))) %>%
    ggplot(aes(fill = PTC, y = n, x = class_code, label = n)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_text(size = 5, position = position_dodge(0.9), color = "black", vjust = -0.2) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_fill_manual(
      values = c("black", "firebrick", "gray50"),
      name = "is_PTC"
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      n.breaks = 5L,
      expand = c(0, 0, 0.1, 0)
    ) +
    facet_wrap(~cds_source) +
    theme_linedraw() +
    theme_Publication() +
    labs(
      title = "NMDtxDB",
      subtitle = "Dataset overview by CDS soruce",
      x = "Reference match",
      y = "Numb. of transcripts"
    )

  ggsave(
    filename = "~/Nextcloud2/nmd_transcriptome_figures/overview_bysource.pdf",
    plot = overview2,
    dpi = "retina"
  )

  ggplot(head(x), aes(fill = PTC, values = n)) +
    geom_waffle(size = .25, n_rows = 10, flip = TRUE) +
    facet_wrap(~class_code, nrow = 1, strip.position = "bottom") +
    scale_x_discrete() +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    ggthemes::scale_fill_tableau(name = NULL) +
    coord_equal() +
    theme_minimal() + #
    theme(panel.grid = element_blank(), axis.ticks.y = element_line()) +
    guides(fill = guide_legend(reverse = TRUE))

  db[["dge"]] %>%
    select(log2FoldChange, padj, contrasts) %>%
    mutate(padj = ifelse(padj <= 1e-20, 1e-20, padj)) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    theme_linedraw() +
    theme_Publication() +
    geom_vline(xintercept = c(-0.1, 0.1), col = "black", linetype = "dotted") +
    geom_hline(yintercept = -log10(1e-20), col = "black", linetype = "dotted") +
    annotate("text", x = 0, y = -log10(1e-20) + 0.5, label = "Max -log10(prob_no_change)") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dotted") +
    labs(
      title = "DGE volcano",
    ) +
    facet_wrap(~contrasts, ncol = 3)
}


plot_volcano <- function(db) {
  tx_anno <- db[["gtf"]] %>%
    filter(type == "transcript", cds_source != "transdecoder") %>%
    select(transcript_id, color) %>%
    left_join(db[["anno"]], by = "transcript_id") %>%
    group_by(gene_id)


  x <- db[["dge"]] %>%
    filter(contrasts == "HEK_NoKO_SMG5KD-KD_Z023-vs-HEK_NoKO_LucKD-KD_Z023") %>%
    left_join(tx_anno, multiple = "all", by = join_by(gene_id, gene_name)) %>%
    distinct() %>%
    mutate(PTC = color == "#FF0000") %>%
    group_by(contrasts, gene_id) %>%
    summarise(log2FoldChange = first(log2FoldChange), PTC = any(PTC), padj = first(padj))

  p2 <- ggplot(x, aes(x = log2FoldChange, fill = PTC)) +
    geom_density(alpha = 0.8) +
    scale_fill_manual(
      values = c("black", "firebrick", "gray50"),
      name = "is_PTC"
    ) +
    xlim(-5, 5) +
    theme_linedraw() +
    theme_Publication() +
    theme(
      plot.margin = margin(0.3, 5.5, 5.5, 5.5)
    )

  p1 <- x %>%
    mutate(padj = ifelse(padj <= 1e-20, 1e-20, padj)) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color = PTC)) +
    geom_point() +
    xlim(-5, 5) +
    scale_color_manual(
      values = c("black", "firebrick", "gray50"),
      name = "is_PTC"
    ) +
    theme_linedraw() +
    theme_Publication() +
    geom_vline(xintercept = c(-0.1, 0.1), col = "black", linetype = "dotted") +
    geom_hline(yintercept = -log10(1e-20), col = "black", linetype = "dotted") +
    annotate("text", x = 0, y = -log10(1e-20) + 0.5, label = "Max -log10(prob_no_change)") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dotted") +
    labs(
      title = "DGE volcano",
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, -1, 5.5)
    )

  p1 / p2 + plot_layout(heights = c(1, .2))
  ggsave(
    filename = "~/Nextcloud2/nmd_transcriptome_figures/volcano.pdf",
    dpi = "retina"
  )
}

plot_summary <- function(variables) {
  ref <- readRDS("phaseFinal/data/gtf.RDS")
  gtf <- rtracklayer::import("phaseFinal/data/gtf_annotated.gtf")
  tx2gene <- readRDS("data/tx2gene.RDS")
  t <- subset(gtf, type == "transcript")

  df <- data.frame(
    name = c(
      "Initial number of transcripts.",
      "Of which are testable (after filtering)",
      "Subset by any CDS overlap",
      "Subset by any CDS overlap (experimental)",
      "Of which have class code not = (novel)",
      "Of which have long-read support",
      "Of which have a PTC",
      "Of which have a PTC and are novel"
    ),
    values = c(
      ref$transcript_id %>% unique() %>% length(),
      nrow(tx2gene %>% filter(keep)),
      t$transcript_id %>% unique() %>% length(),
      subset(t, cds_source != "transdecoder")$transcript_id %>% unique() %>% length(),
      subset(t, cds_source != "transdecoder" & class_code != "=")$transcript_id %>% unique() %>% length(),
      subset(t, cds_source != "transdecoder" & lr_support == "TRUE")$transcript_id %>% unique() %>% length(),
      subset(t, cds_source != "transdecoder" & color == "#FF0000")$transcript_id %>% unique() %>% length(),
      subset(t, cds_source != "transdecoder" & color == "#FF0000" & class_code != "=")$transcript_id %>% unique() %>% length()
    )
  )

  df <- df %>% mutate(name = fct_reorder(name, values, .desc = FALSE))

  ggplot(df, aes(name, y = values)) +
    geom_bar(stat = "identity") +
    theme_light(16) +
    labs(x = element_blank(), y = element_blank()) +
    scale_y_log10(
      "",
      limits = c(1, 1e6),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    coord_flip() +
    geom_text(aes(label = scales::comma(values)), vjust = 0.5, hjust = -0.1)
  ggsave("phaseFinal/nmd_transcriptome_stats.pdf", width = 10, height = 7, dpi = 300)

  x <- mcols(t) %>%
    as_tibble() %>%
    mutate(ptc = ifelse(color == "#FF0000", T, F)) %>%
    select(transcript_id, ptc, cds_source, lr_support, class_code) %>%
    mutate(
      cds_source = factor(
        cds_source,
        levels = c("ensembl", "hek293gao", "ribotish", "openprot", "transdecoder")
      ),
      class_code = factor(
        if_else(class_code != "=", "novel", "not_novel"),
        levels = c("novel", "not_novel")
      ),
      ptc = factor(ptc, levels = c("TRUE", "FALSE")),
      lr_support = factor(lr_support, levels = c("TRUE", "FALSE")),
      transcript_id = str_replace(transcript_id, "_.*", "")
    ) %>%
    arrange(ptc, cds_source, class_code, lr_support) %>%
    filter(cds_source != "transdecoder")

  ggsave(
    filename = "~/Nextcloud2/nmd_transcriptome_figures/volcano.pdf",
    dpi = "retina"
  )


  # x2 <- x %>%
  #   sample_frac(.10) %>%
  #   group_by(transcript_id) %>%
  #   summarise(
  #     across(everything(), first))  %>%
  #   ungroup() %>%
  #   group_by(across(-transcript_id)) %>%
  #   tally()
  #
  # x2 <- x2 %>%
  #   ungroup() %>%
  #   mutate(across(where(is.factor), as.character))
  #
  # PieDonut(x2, aes(ptc, cds_source, class_code, lr_support, count=n), title = "nmdTXdb stats")
}


plot_tc_dist_ref <- function() {
  # suppressPackageStartupMessages({
  #   library(dplyr)
  #   library(plyranges)
  #   library(ORFik)
  #   library(rtracklayer)
  #   library(Biostrings)
  #   library(GenomicFeatures)
  # })
  #
  # ref <- import("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.SIRV.gtf")
  #
  # ungap <- function(x) {
  #   stopifnot(inherits(x, "IRanges"))
  #   IRanges(end=cumsum(width(x)), width = width(x))
  # }
  #
  # cds <- ref %>%
  #   filter(type == "CDS") %>%
  #   plyranges::select(transcript_id)
  #
  # gr <- ref %>%  # change to tx
  #   filter(type == 'exon') %>%
  #   plyranges::select(transcript_id) %>%
  #   filter(transcript_id %in% unique(cds$transcript_id) )
  #
  # grl <- gr %>% split(., ~transcript_id) # tx_l
  # x_range <- range(grl)
  # if (any(elementNROWS(x_range) != 1L))
  #   stop("Empty or multi-strand/seqname elements not supported by BED")
  #
  # cds <- pmapToTranscripts(cds, grl[cds$transcript_id])
  # gr <- pmapToTranscripts(gr, grl[gr$transcript_id])
  # cds <- split(ranges(cds), seqnames(cds))
  # gr <- split(ranges(gr), seqnames(gr))
  # stopifnot(all(names(gr) == names(cds)))
  #
  # # five prime plus cds
  # fcds <- cds
  # start(fcds) <- gr %>% heads(1) %>% start() # here we extend the start to 1
  # fcds <- range(fcds)
  # fcdsl <- intersect(fcds, gr)
  #
  # gr_ungap <- endoapply(gr, ungap)
  # fcds_ungap <- endoapply(fcdsl, ungap)
  # stop <- fcds_ungap %>% tails(1) %>% end()
  # last_ejc <- gr_ungap %>% tails(2) %>% heads(1) %>% end()
  # dist <- last_ejc - 1 - stop
  #
  # library(ggplot2)
  #
  # df <- ref %>%
  #   filter(type == 'transcript') %>%
  #   as.data.frame() %>%
  #   dplyr::select(transcript_id, strand, transcript_biotype) %>%
  #   filter(transcript_biotype %in% c('nonsense_mediated_decay', 'protein_coding')) %>%
  #   mutate(dist = unlist(dist, use.names = FALSE)[match(transcript_id, names(dist))])
  #
  # df %>%
  #   ggplot() +
  #   geom_density(aes(x = dist, fill=transcript_biotype), alpha=.4, color='white') +
  #   geom_vline(xintercept=50, linetype='dashed') +
  #   facet_wrap(~strand) +
  #   coord_cartesian(xlim=c(-1000, 1000))
  # ggsave('cds_task_refonly.pdf', height = 7, width = 10)
}
