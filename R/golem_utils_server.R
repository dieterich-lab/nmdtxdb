#' Inverted versions of in, is.null and is.na
#'
#' @noRd
#'
#' @examples
#' 1 %not_in% 1:10
#' not_null(NULL)
`%not_in%` <- Negate(`%in%`)

not_null <- Negate(is.null)

not_na <- Negate(is.na)

#' Removes the null from a vector
#'
#' @noRd
#'
#' @example
#' drop_nulls(list(1, NULL, 2))
drop_nulls <- function(x) {
  x[!sapply(x, is.null)]
}

#' If x is `NULL`, return y, otherwise return x
#'
#' @param x,y Two elements to test, one potentially `NULL`
#'
#' @noRd
#'
#' @examples
#' NULL %||% 1
"%||%" <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

#' If x is `NA`, return y, otherwise return x
#'
#' @param x,y Two elements to test, one potentially `NA`
#'
#' @noRd
#'
#' @examples
#' NA %|NA|% 1
"%|NA|%" <- function(x, y) {
  if (is.na(x)) {
    y
  } else {
    x
  }
}

#' Typing reactiveValues is too long
#'
#' @inheritParams reactiveValues
#' @inheritParams reactiveValuesToList
#'
#' @noRd
rv <- function(...) shiny::reactiveValues(...)
rvtl <- function(...) shiny::reactiveValuesToList(...)




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


log10_or_max <- function(x) {
  log10x <- log10(x)
  log10x[log10x <= -Inf] <- min(log10x[is.finite(log10x)])
  log10x[log10x >= Inf] <- max(log10x[is.finite(log10x)])
  log10x
}

#' @importFrom httr modify_url GET http_type http_error status_code accept content
#' @importFrom jsonlite fromJSON
get_gene_info <- function(gene_id) {
  url <- modify_url(
    url = "https://mygene.info",
    path = paste0("v3/gene/", gene_id),
    query = "fields=symbol,summary,alias,uniprot,name"
  )

  resp <- GET(url, accept("application/json"))

  if (http_type(resp) != "application/json") {
    stop("API did not return json", call. = FALSE)
  }

  parsed <- fromJSON(content(resp, "text"), simplifyVector = TRUE)
  if (http_error(resp)) {
    stop(
      sprintf(
        "mygene.info API request failed [%s]\n%s\n<%s>",
        status_code(resp),
        parsed$message,
        parsed$documentation_url
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      content = parsed,
      url = url,
      response = resp
    ),
    class = "mygene_info_api"
  )
}

#' @importFrom htmltools tagAppendAttributes
#' @importFrom stringr str_interp str_glue
#' @note remove db calls and test
render_gene_card <- function(gene_id) {
  parsed <- get_gene_info(gene_id)$content

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
          str_interp("Also known as: ${paste(parsed$alias, collapse=', ')}")
        )
      ),
      div(
        "Ensembl:",
        a(
          gene_id,
          href = paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", gene_id),
          target = "_blank"
        )
      ),
      conditionalPanel(
        !is.null(parsed$uniprot$`Swiss-Prot`),
        p(
          "Uniprot:",
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

#' Plot a transcript on the cDNA coordiantes highlighting PTC containing
#' transcripts
#' @import dplyr ggtranscript
#'
plot_annotation_cdna <- function(bed12) {
  feat_colors <- c("TRUE" = "firebrick", "FALSE" = "black")
  exon <- bed12$cdna_blocks
  names(exon) <- bed12$cdna_thick$name
  exon <- exon %>% dplyr::bind_rows(.id = "cds_id")

  cds <- bed12$cdna_thick %>%
    dplyr::bind_rows() %>%
    group_by(names) %>%
    summarize(start = min(start), end = max(end)) %>%
    ungroup() %>%
    dplyr::rename(cds_id=names) %>%
    left_join(bed12 %>% select(cds_id, is_ptc), by = 'cds_id')

  text <- exon %>%
    group_by(cds_id) %>%
    mutate(eid = 1: n(), x = (start + end) / 2)

  exon %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = cds_id
    )) +
    geom_range(
      fill = "white",
      height = 0.25
    ) +
    geom_range(
      data = cds,
      height = 0.40,
      alpha = .50,
      aes(
        fill = is_ptc
      )
    ) +
    scale_fill_manual(values = feat_colors) +
    geom_text(
      data = text,
      size = 3,
      vjust = 1.5,
      aes(label = eid, x = x)
    ) +
    labs(fill = "PTC:") +
    theme_linedraw() +
    labs(y = "", x = "")
}


#' Simplifies grid creation
#' @note does not support mobile grid.
#' @importFrom shiny.semantic grid_template
create_grid <- function(areas, cols_width, rows_height) {
  shiny.semantic::grid_template(
    default = list(
      areas = areas,
      cols_width = cols_width,
      rows_height = rows_height
    )
  )
}

#' Extract Genomic Position from GTF
#'
#' This function extracts the genomic position from a GTF file by grouping the data based on a specified column and summarizing the sequence name, start, and end coordinates.
#'
#' @param gtf A data frame containing the GTF data.
#' @param col The column to group by for summarization.
#'
#' @return A character string representing the genomic position in the format "chr{seqnames}:{start}-{end}".
#'
#' @examples
#' gtf <- data.frame(
#'   seqnames = c("chr1", "chr1", "chr2", "chr2"),
#'   start = c(100, 200, 300, 400),
#'   end = c(500, 600, 700, 800),
#'   gene_id = c("Gene1", "Gene2", "Gene3", "Gene4")
#' )
#' position_from_gtf(gtf, "gene_id")
#'
#' @importFrom dplyr group_by summarise first
#' @importFrom stringr str_glue_data
#'
position_from_gtf <- function(gtf) {
  gtf %>%
    mutate(chr = ifelse(str_detect(seqnames, "chr"), "", "chr")) %>%
    str_glue_data("{chr}{seqnames}:{start}-{end}")
}


#' Make Unique Character Vectors with Suffixes
#'
#' Takes a character vector and makes each element unique by appending an increment
#'
#' @param strings A character vector with potentially non-unique elements.
#' @return A character vector where each element has been made unique by appending
#' an underscore and an occurrence number.
#' @examples
#' original_vector <- c("apple", "banana", "apple", "orange", "banana", "banana")
#' make_unique(original_vector)
#' @export
make_unique <- function(strings) {
  counts <- ave(strings, strings, FUN = seq_along)
  return(paste0(strings, "_", counts))
}


#' Create Trackhub URL
#'
#' This function generates a URL for accessing a track hub on the UCSC Genome Browser.
#'
#' @param base_url The base URL of the UCSC Genome Browser. Default is "http://genome-euro.ucsc.edu/".
#' @param db The genome database. Default is "hg38".
#' @param position The genomic position. Default is NA for no posiitoj.
#' @param hub The URL of the track hub. Default is "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/nmd_transcriptome.hub.txt".
#'
#' @return A URL for accessing the specified track hub on the UCSC Genome Browser.
#'
#' @examples
#' create_trackhub_url()
#' create_trackhub_url(db = "mm10", position = "chr1:1000-2000")
#' create_trackhub_url(db = "mm10", position = "chr1:1000-2000", hub = NA)
#' @importFrom httr modify_url
#'
create_trackhub_url <- function(base_url = "http://genome-euro.ucsc.edu/", db = "hg38", position = NA, hub = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/nmd_transcriptome.hub.txt") {
  query <- list(
    db = db,
    position = position,
    hubUrl = hub
  )

  modify_url(
    url = base_url,
    path = "cgi-bin/hgTracks",
    query = query[!is.na(query)]
  )
}

load_metadata <- function(db) {
  db[["metadata"]] %>%
    select(contrasts, Knockdown, Knockout, cellline, clone) %>%
    filter(Knockdown != "LucKD") %>%
    distinct() %>%
    mutate(
      name = str_replace(contrasts, "-vs-.*", ""),
      Knockout = str_replace(Knockout, "_", ""),
      clone = str_replace_na(""),
      contrast_label = str_glue_data(
        .,
        "<b> {Knockdown} </b> <br> <i>Cell-line</i>: {cellline}; <i>KO</i>: {Knockout} <i>clone</i>: {clone}"
      ),
      label = str_c(Knockdown, Knockout, cellline, clone,  sep = "_")
    ) %>% View()
}

with_tooltip <- function(value, tooltip) {
  tags$abbr(
    style = "text-decoration: underline; text-decoration-style: dotted; cursor: help",
    title = tooltip, value
  )
}


#' Plot Boxplot with Gene Highlight
#'
#' This function generates a boxplot with a highlighted gene based on the specified contrast.
#'
#' @param data A data frame containing the necessary columns: gene_name, contrasts, log2FoldChange.
#' @param gene_name The name of the gene to highlight.
#' @param xlims minimum, maximum values for x-axis
#'
#'
#' @return A ggplot object displaying the boxplot with the highlighted gene.
#'
#' @examples
# data <- data.frame(gene_name = c("Gene1", "Gene2", "Gene2"),
#                    contrasts = c("Contrast1", "Contrast2", "Contrast2"),
#                    log2FoldChange = c(1.5, -0.8, 2.2))
# fc_boxplot(data, "Gene2", "Contrast2")
#'
#' @import dplyr
#' @import ggplot2
fc_boxplot <- function(data, index, gene_l2fc, xlims) {
  data <- data[index, ]
  gene_name <- data[[1, "gene_name"]]
  contrast <- data[[1, "contrasts"]]
  data <- mutate(data, y = 1, label = gene_name)

  gene_l2fc <- gene_l2fc %>% filter(contrasts == !!contrast)
  text_color <- ifelse(data$log2FoldChange > 0, "red", "blue")
  p <- gene_l2fc %>%
    ggplot(aes(y = contrasts, x = log2FoldChange)) +
    geom_boxplot() +
    geom_text(
      data = data,
      aes(label = label, y = y),
      position = position_nudge(y = 0.08),
      hjust = 0,
      colour = text_color,
      size = 3.5
    ) +
    geom_vline(data = data, colour = text_color, aes(xintercept = log2FoldChange)) +
    lims(x = xlims) +
    labs(y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank()
    )
  return(p)
}
