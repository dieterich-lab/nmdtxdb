grid <- create_grid(
  rbind(
    c("top"),
    c("bottom_left")
  ),
  c("100%", "100%"),
  c("70%", "30%")
)


#' transcript view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS
#' @importFrom reactable reactableOutput
#' @importFrom patchwork plot_layout
mod_transcript_ui <- function(id) {
  ns <- NS(id)

  grid(
    grid,
    top = reactableOutput(ns("table_transcript")) %>%
      shinycssloaders::withSpinner(),
    bottom_left = plotOutput(ns("gene_structure")) %>% shinycssloaders::withSpinner(),
    bottom_right = NULL
  )
}


#' Custom GGPlot Function
#'
#' This function takes a dataframe and a vector of \code{cds_id} values, and creates a customized GGPlot using ggplot2.
#'
#' @param df A dataframe containing the necessary columns.
#' @param y_labs A vector of \code{cds_id} values to filter the dataframe.
#'
#' @return A GGPlot object.
#'
#' @import ggplot2
#' @importFrom dplyr select filter mutate
#' @importFrom ggplot2 labs geom_point theme_linedraw facet_grid scale_color_gradient2
#' @importFrom ggplot2 element_text element_blank
#'
#' @examples
# df <- data.frame(
#   cds_id = c("id1", "id2", "id3"),
#   cds_source2 = c("source1", "source2", "source1"),
#   label = c("label1", "label2", "label3"),
#   padj = c(0.05, 0.01, 0.001),
#   log2fold = c(2, -1, 3)
# )
# y_labs <- c("id1", "id3")
# build_dotplot(df, y_labs)
#'
build_dotplot <- function(df, y_labs) {
  df %>%
    select(cds_id, cds_source2, label, padj, log2fold) %>%
    filter(cds_id %in% y_labs) %>%
    mutate(cds_id = factor(cds_id, levels = y_labs), label = str_replace_all(label, "_", "\n")) %>%
    ggplot(aes(x = label, y = cds_id, size = -log10(as.numeric(padj)), color = log2fold)) +
    labs(color = "log2fold: ", size = "-log10(padj):") +
    geom_point() +
    theme_linedraw() +
    labs(y = "", x = "") +
    facet_grid(cds_source2 ~ label, scales = "free") +
    scale_color_gradient2(low = "black", mid = "aliceblue", high = "firebrick") +
    theme(
      text = element_text(size = 10),
      axis.text.y = element_blank(),
      axis.ticks.length.y = unit(0, "pt"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_text(size = 5),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
}


#' Create Transcript Plot
#'
#' This function takes a dataframe containing GTF information and creates a
#' customized transcript plot using ggplot2.
#'
#' @param gtf A dataframe containing GTF information.
#'
#' @return A GGPlot object.
#'
#' @import ggplot2
#' @importFrom ggplot2 labs scale_fill_manual facet_wrap theme_linedraw
#' @importFrom ggtranscript geom_range geom_intron
#' @importFrom grid arrow unit
#' @importFrom ggplot2 element_text element_blank
#'
#' @examples
#' gtf <- data.frame(
#'   type = c("exon", "CDS", "exon", "exon", "CDS", "exon"),
#'   start = c(100, 150, 200, 250, 300, 350),
#'   end = c(130, 170, 220, 260, 310, 360),
#'   Name = c("transcript1", "transcript1", "transcript2", "transcript3", "transcript2", "transcript4"),
#'   cds_source2 = c("source1", "source1", "source2", "source2", "source3", "source3"),
#'   PTC = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
#'   strand = c("+", "-", "+", "+", "-", "+")
#' )
#' build_transcript_plot(gtf)
#'
build_transcript_plot <- function(gtf) {
  if (nrow(gtf) < 1) {
    return("No CDS source to show.")
  }

  exons <- gtf %>% filter(type == "exon")
  cds <- gtf %>% filter(type == "CDS")
  introns <- to_intron(exons, group_var = "Name")
  feat_colors <- c("TRUE" = "firebrick", "FALSE" = "black")

  exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = Name
    )) +
    geom_range(
      aes(
        fill = PTC,
        height = 0.25
      )
    ) +
    geom_range(
      data = cds,
      aes(
        fill = PTC
      )
    ) +
    geom_intron(
      data = introns,
      aes(strand = strand),
      arrow = grid::arrow(ends = "last", length = grid::unit(0.4, "lines"))
    ) +
    scale_fill_manual(values = feat_colors) +
    labs(fill = "PTC:") +
    facet_wrap(~cds_source2, ncol = 1, scales = "free_y", strip.position = "left") +
    theme_linedraw() +
    labs(y = "") +
    theme(
      text = element_text(size = 10),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "top"
    )
}

#' transcript view Server Functions
#' @import dplyr
#' @importFrom purrr map
#' @importFrom reactable getReactableState renderReactable
#' @import stringr
#' @importFrom scales scientific
#' @noRd
mod_transcript_server <- function(id, conn, tx, contrast, cds) {
  moduleServer(id, function(input, output, session) {
    l2fc <- tbl(conn, "dte") %>%
      filter(contrasts %in% !!contrast) %>%
      collect()

    transcripts <- tbl(conn, "gtf") %>%
      filter(transcript_id %in% !!tx, type == "transcript", cds_source %in% !!cds) %>%
      select(Name, transcript_id, cds_source, color, seqnames, start, end) %>%
      collect() %>%
      mutate(
        PTC = color == "#FF0000",
        position = position_from_gtf(.),
        trackhub_url = purrr::map(
          position,
          \(x)  create_trackhub_url(position = x)
        ) %>% unlist()
      )

    not_transcrips <- tbl(conn, "gtf") %>%
      filter(transcript_id %in% !!tx, type != "transcript") %>%
      collect() %>%
      left_join(
        x = transcripts %>% select(Name, cds_source, PTC),
        y = select(., !c(cds_source)),
        by = "Name", multiple = "all"
      ) %>%
      left_join(cds_source_choices, by = "cds_source", multiple = "all")

    dte <- tbl(conn, "dte") %>%
      filter(
        transcript_id %in% !!tx,
        contrasts %in% !!contrast
      ) %>%
      select(transcript_id, contrasts, padj, log2fold) %>%
      collect() %>%
      left_join(load_metadata(conn), by = "contrasts") %>%
      select(-name) %>%
      collect()

    anno <- tbl(conn, "anno") %>%
      filter(transcript_id %in% !!tx) %>%
      select(
        ref_gene_name, transcript_id, ref_transcript_name,
        ref_transcript_id, lr_support, class_code
      ) %>%
      collect()

    cds_position <- not_transcrips %>%
      filter(type == "CDS") %>%
      group_by(transcript_id, cds_id) %>%
      summarise(seqnames = first(seqnames), start = min(start), end = max(end)) %>%
      ungroup() %>%
      mutate(cds_position = position_from_gtf(.)) %>%
      select(transcript_id, cds_id, cds_position) %>%
      mutate(
        cds_trackhub_url = purrr::map(
          cds_position,
          \(x) create_trackhub_url(position = x)
        ) %>% unlist()
      )

    cds_data <- transcripts %>%
      rename(cds_id = Name) %>%
      left_join(cds_position, by = c("transcript_id", "cds_id")) %>%
      mutate(
        cds_source2 = case_when(
          cds_source == "ensembl" ~ "Ensembl",
          cds_source == "hek293gao" ~ "Gao et al., 2015",
          cds_source == "openprot" ~ "OpenProt",
          cds_source == "ribotish" ~ "Zhang et al., 2017",
          cds_source == "transdecoder" ~ "TransDecoder"
        )
      )

    df <- anno %>%
      left_join(dte, by = "transcript_id", multiple = "all") %>%
      left_join(cds_data, by = "transcript_id", multiple = "all") %>%
      select(transcript_id, ref_transcript_name, PTC, lr_support, everything()) %>%
      mutate(
        log2fold = round(log2fold, 2),
        padj = padj %>% scientific(),
        class_code = case_when(
          class_code == "=" ~ "Complete",
          class_code == "n" ~ "Retained introns",
          class_code %in% c("c", "j", "k") ~ "Splicing variants",
          TRUE ~ "Other"
        ),
      ) %>%
      select(-c(seqnames, start, end, color)) %>%
      distinct(transcript_id, cds_id, cds_source, contrasts, .keep_all = TRUE)


    output$table_transcript <- renderReactable({

      reactable(
        df,
        defaultSorted = c("PTC"),
        defaultPageSize = 5,
        bordered = TRUE,
        highlight = TRUE,
        wrap = FALSE,
        details = function(index) {
          if (!is.na(df[[index, "contrasts"]])) {
            data <- df[index, ]
            text_color <- ifelse(data$log2fold > 0, "red", "blue")

            p <- l2fc %>%
              filter(!is.na(log2fold), contrasts == data[["contrasts"]]) %>%
              ggplot(aes(y = contrasts, x = log2fold)) +
              geom_boxplot() +
              geom_vline(data = data, colour = text_color, aes(xintercept = log2fold)) +
              labs(y = "") +
              theme_minimal() +
              theme(
                axis.text.y = element_blank()
              )

            htmltools::plotTag(p, alt = "plots", height = 150)
          }
        },
        theme = reactableTheme(
          borderColor = "#dfe2e5",
          stripedColor = "#f6f8fa",
          highlightColor = "#FFFFBF",
          cellPadding = "8px 12px"
        ),
        defaultColDef = colDef(
          sortNALast = TRUE
        ),
        columns = list(
          cds_id = colDef(
            header = with_tooltip(
              "CDS_info", "CDS ID, provenance and URL to the trackhub."
            ),
            width = 210,
            show = TRUE,
            align = "left",
            vAlign = "center",
            html = TRUE,
            cell = JS("
function (cellInfo) {
  const cds_source = cellInfo.row['cds_source2'];
  const cds_id = cellInfo.row['cds_id'];
  const cds_trackhub_url = cellInfo.row['cds_trackhub_url'];

  if (cds_id === null) {
    return 'No CDS';
  } else {
    return `<div><a href='${cds_trackhub_url}' target='_blank'>${cds_id}</a></div><div><small><i>Source</i>: ${cds_source}</small></div>`;
  }
}")
          ),
          PTC = colDef(
            header = with_tooltip(
              "PTC", "\u2713 if the transcript has a PTC else \u274c."
            ),
            defaultSortOrder = "desc",
            show = TRUE,
            width = 80,
            align = "center",
            vAlign = "center",
            cell = function(value) {
              if_else(value == FALSE, "\u274c", "\u2713", missing = "")
            }
          ),
          lr_support = colDef(
            header = with_tooltip(
              "LRS", "\u2713 if the transcript has long read support else \u274c."
            ),
            width = 60,
            show = TRUE,
            align = "center",
            vAlign = "center",
            cell = function(value) {
              if_else(value == "FALSE", "\u274c", "\u2713", missing = "")
            }
          ),
          transcript_id = colDef(
            header = with_tooltip(
              "transcript_id", "Identifier and link to the UCSC Genome Browser Trackhub."
            ),
            width = 160,
            show = TRUE,
            vAlign = "center",
            html = TRUE,
            cell = JS("
function (cellInfo) {
  const tid = cellInfo.row['transcript_id'];
  const url = cellInfo.row['trackhub_url'] || undefined;

  if (url === undefined) {
    return tid;
  } else {
    return `<a href='${url}' target='_blank'>${tid}</a>`;
  }
}")
          ),
          trackhub_url = colDef(
            show = FALSE,
          ),
          ref_transcript_name = colDef(
            header = with_tooltip(
              "ref_transcript", "Ensembl transcript name followed by type of match between reference and assembly."
            ),
            width = 160,
            vAlign = "center",
            show = TRUE,
            html = TRUE,
            cell = JS("
function(cellInfo) {
  const cc = cellInfo.row['class_code'];
  const ti = cellInfo.row['ref_transcript_id'];
  const lab1 = '<div>' + '<strong>' + cellInfo.value + '</strong> </div>'
  const lab2 = '<div><small><i>Match</i>: ' + cc + '</small></div>';

  const url = 'http://www.ensembl.org/id/' + ti;
  return `<a href='${url}' target='_blank'>${lab1}</a>${lab2}`;

}")
          ),
          contrasts = colDef(
            header = with_tooltip(
              "contrasts", "DTE comparison in the format treatment vs control."
            ),
            width = 160,
            show = TRUE,
            html = TRUE,
            vAlign = "center",
            cell = JS("
function(cellInfo) {
  const kd = cellInfo.row['Knockdown'];
  const cl = cellInfo.row['cellline'];
  const ko = cellInfo.row['Knockout'] || 'NoKO';

  if (kd === null) {
    return 'Not tested';
  } else {
    return `<div><strong>${kd}</strong><br><small><i>Cell-line</i>: ${cl};<i> KO</i>: ${ko}</small></div>`;
  }
}")
          ),
          padj = colDef(
            header = with_tooltip(
              "padj", "DTE adjusted p-value."
            ),
            filterable = FALSE,
            show = TRUE,
            width = 100,
            align = "right"
          ),
          log2fold = colDef(
            header = with_tooltip(
              "log2fc", "DTE log2FoldChange."
            ),
            show = TRUE,
            width = 80,
            format = colFormat(digits = 2),
            filterable = FALSE
          ),
          Knockdown = colDef(
            show = FALSE
          ),
          cds_source = colDef(
            show = FALSE
          ),
          cds_source2 = colDef(
            show = FALSE
          ),
          cds_id = colDef(
            show = FALSE
          ),
          cds_position = colDef(
            show = FALSE
          ),
          cds_trackhub_url = colDef(
            show = FALSE
          ),
          Knockout = colDef(
            show = FALSE
          ),
          cellline = colDef(
            show = FALSE
          ),
          contrast_label = colDef(
            show = FALSE
          ),
          label = colDef(
            show = FALSE
          ),
          class_code = colDef(
            show = FALSE
          ),
          ref_gene_name = colDef(
            show = FALSE
          ),
          ref_transcript_name = colDef(
            show = FALSE
          ),
          ref_transcript_id = colDef(
            show = FALSE
          ),
          position = colDef(
            show = FALSE
          )
        )
      )
    })
    output$gene_structure <- renderPlot(
      {
        p1 <- build_transcript_plot(not_transcrips)
        y_labs <- lapply(
          ggplot_build(p1)$layout$panel_params,
          \(x) x$y$get_labels()) |>
          unlist()
        p2 <- build_dotplot(df, y_labs)
        p1 + p2 +
          patchwork::plot_layout(widths = c(4, 1), guides = "collec") &
          theme(
            legend.position = "top",
            legend.text = element_text(size=8),
            legend.key.height = unit(0.1, "cm"),
            panel.spacing.x= unit(0.1, "cm"),
            panel.spacing.y= unit(0.1, "cm")
          )
      },
      res = 150
    )
  })
}
