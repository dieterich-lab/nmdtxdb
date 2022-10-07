
#' @importFrom shiny column
col_6 <- function(...) {
  column(6, ...)
}

#' @importFrom shiny column
col_3 <- function(...) {
  column(3, ...)
}

#' @importFrom shiny column
col_2 <- function(...) {
  column(2, ...)
}

#' @importFrom shiny column
col_5 <- function(...) {
  column(5, ...)
}


adv_grid <- grid_template(
  default = list(
    areas = rbind(
      c("top_left", "top_left", "top_right"),
      c("bottom_left", "bottom_mid", "bottom_right")
    ),
    rows_height = c("30%", "70%"),
    cols_width = c("20%", "30%", "50%")
  )
)


bigwig <- list(
  control1 = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/hg38/33F-1.bigWig",
  SMG6kd_SMG7ko = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/hg38/33F-14.bigWig"
)

log10_or_max <- function(x) {
  log10x <- log10(x)
  log10x[log10x <= -Inf] <- min(log10x[is.finite(log10x)])
  log10x[log10x >= Inf] <- max(log10x[is.finite(log10x)])
  log10x
}

get_gene_info <- function(gene_id) {
  url <- httr::modify_url(
    url = "https://mygene.info",
    path = paste0("v3/gene/", gene_id),
    query = "fields=symbol,summary,alias,uniprot,name"
  )

  resp <- httr::GET(url, httr::accept("application/json"))
  jsonlite::fromJSON(httr::content(resp, "text"), simplifyVector = TRUE)
}


render_gene_card <- function(gene_id) {
  parsed <- get_gene_info(gene_id)
  transcripts <- tbl(conn, "gtf") %>%
    filter(type == "transcript", gene_id == !!gene_id) %>%
    collect()

  dge_l2fc <- tbl(conn, "dge") %>%
    filter(gene_id == !!gene_id) %>%
    collect() %>%
    pull(log2FoldChange)
  is_nmd <- transcripts %>% with(table(transcript_biotype == "nonsense_mediated_decay"))
  n_transcripts <- nrow(transcripts)
  n_nmd <- ifelse(length(is_nmd) == 2, is_nmd[2], 0)
  has_support <- tbl(conn, "has_support") %>% collect()
  n_support <- length(intersect(transcripts$transcript_id, has_support$transcript_id))

  card(
    align = "left",
    div(
      a(class = "ui blue ribbon label", "Gene information"),
      class = "content",
      br(),
      br(),
      div(class = "header", parsed$symbol),
      div(parsed$name),
      hr(),
      div(
        class = "meta",
        HTML(
          stringr::str_interp("Also known as: ${paste(parsed$alias, collapse=', ')}")
        )
      ),
      div(str_glue("NMD isoforms: {n_nmd}/{sum(is_nmd)}")),
      div(
        str_glue("Up-regulated in {sum(dge_l2fc > 0.1, na.rm=TRUE)}/{n_transcripts} cohorts"),
        icon("question circle"),
      ) %>% htmltools::tagAppendAttributes(., "data-tooltip" = "Number of datasets with l2fc > 0.1"),
      # p("Novel transcripts : X/XX"),
      div(
        str_glue("Long read evidence: {n_support}/{n_transcripts}"),
        icon("question circle"),
      ) %>% htmltools::tagAppendAttributes(
        ., "data-tooltip" = "with identical intron chain"),
      div(
        "Ensembl: ",
        a(
          gene_id,
          href = paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", gene_id),
          target = "_blank"
        )
      ),
      conditionalPanel(
        !is.null(parsed$uniprot$`Swiss-Prot`),
        p(
          "Uniprot: ",
          a(
            parsed$uniprot$`Swiss-Prot`,
            href = paste0("https://www.uniprot.org/uniprotkb/", parsed$uniprot$`Swiss-Prot`),
            target = "_blank"
          )
        )
      )
    )
  )
}


send_toast <- function(msg, session, position = "top right", class = "warning", icon = "exclamation") {
  toast(
    msg,
    class = class,
    session = session,
    toast_tags = list(position = position, showIcon = icon)
  )
}

plot_coverage <- function(gtf) {
  region <- GenomicRanges::reduce(GenomicRanges::GRanges(gtf))
  Signac::BigwigTrack(region, bigwig) +
    theme_minimal()
}

plot_annotation <- function(gtf) {
  gtf <- gtf %>%
    mutate(name = paste0(transcript_name, " ", transcript_biotype))

  exons <- gtf %>% dplyr::filter(type == "exon")
  cds <- gtf %>% dplyr::filter(type == "CDS")
  introns <- ggtranscript::to_intron(exons, group_var = "transcript_id")

  p1 <- exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = name
    )) +
    ggtranscript::geom_range(
      aes(
        height = 0.25
      )
    ) +
    ggtranscript::geom_range(
      data = cds
    ) +
    ggtranscript::geom_intron(
      data = introns,
      aes(strand = strand)
    ) +
    theme_minimal() +
    labs(y = "") +
    theme(
      axis.ticks = element_blank(),
      legend.position = c(0.87, 0.75)
    )
}

gene_view_grid <- grid_template(
  default = list(
    areas = rbind(c("left", "right")),
    cols_width = c("30%", "auto"),
    rows_height = c("auto")
  )
)

conn <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "nmd_transcriptome",
  host = "***REMOVED***",
  port = ***REMOVED***,
  password = "***REMOVED***",
  user = "***REMOVED***"
)
