this_grid <- create_grid(
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
#' @importFrom shiny.semantic grid
#' @importFrom shinycssloaders withSpinner
#' @importFrom reactable reactableOutput
mod_transcript_ui <- function(id) {
  ns <- NS(id)

  grid(
    this_grid,
    top = reactableOutput(ns("table_transcript")) %>%
      withSpinner(),
    bottom_left = div(
      plotOutput(ns("gene_structure")) %>%
        withSpinner(),
      downloadButton(ns("downloadPlot"), "Download Plot")
    )
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
#' df <- data.frame(
#'   cds_id = c("id1", "id2", "id3"),
#'   source = c("source1", "source2", "source1"),
#'   label = c("label1", "label2", "label3"),
#'   padj = c(0.05, 0.01, 0.001),
#'   log2fold = c(2, -1, 3)
#' )
#' y_labs <- c("id1", "id3")
#' build_dotplot(df, y_labs)
#'
build_dotplot <- function(df, y_labs) {
  df <- df %>%
    select(cds_id, name, padj, log2fold) %>%
    filter(cds_id %in% y_labs) %>%
    mutate(label = str_replace_all(name, '_', '\n'))

  ggplot(
    df,
    aes(
      x = label,
      y = cds_id,
      size = -log10(as.numeric(padj) + 1e-20),
      color = log2fold
    )
  ) +
    labs(color = "log2fold: ", size = "-log10(padj): ") +
    geom_point() +
    theme_linedraw() +
    labs(y = "", x = "") +
    scale_size_continuous(limits = c(0, 20), range = c(0, 6)) +
    scale_x_discrete(
      limits=na.omit(unique(df$label)),
      guide = guide_axis(n.dodge=2)) +
    scale_color_gradient2(
      low = "blue", mid = "lightyellow", high = "red",
      oob = scales::squish, limits = c(-5, 5)
    ) +
    theme(
      text = element_text(size = 12),
      axis.text.y = element_blank(),
      axis.ticks.length.y = unit(0, "pt"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 5),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
}


#' transcript view Server Functions
#' @import dplyr
#' @importFrom purrr map
#' @importFrom reactable getReactableState renderReactable
#' @importFrom patchwork plot_layout
#' @import stringr
#' @importFrom scales scientific
#' @noRd
mod_transcript_server <- function(id, db, tx, contrast, cds) {
  moduleServer(id, function(input, output, session) {
    l2fc <- db[["dte"]] %>%
      filter(contrasts %in% !!contrast)

    transcripts <- db[["bed12"]] %>%
      filter(name %in% !!tx, source %in% !!cds) %>%
      mutate(
        cds_id = .$cdna_thick$name,
        position = glue::glue("chr{seqnames}:{start}-{end}"),
        cds_position = glue::glue("chr{seqnames}:{thick}"),
        trackhub_url = purrr::map(
          position,
          \(x) create_trackhub_url(position = x)
        ) %>% unlist(),
        cds_trackhub_url = purrr::map(
          cds_position,
          \(x) create_trackhub_url(position = x)
        ) %>% unlist()
      )

    dte <- db[["dte"]] %>%
      filter(
        transcript_id %in% !!tx,
        contrasts %in% !!contrast
      ) %>%
      select(transcript_id, contrasts, padj, log2fold) %>%
      left_join(load_metadata(db), by = "contrasts")

    anno <- db[["anno"]] %>%
      filter(transcript_id %in% !!tx) %>%
      select(
        ref_gene_name, transcript_id, ref_transcript_name,
        ref_transcript_id, lr_support, match
      )

    fname <- file.path(tempdir(), paste0("nmdtxdb_", anno$ref_gene_name[1], ".png"))

    df <- anno %>%
      left_join(dte, by = "transcript_id", multiple = "all") %>%
      left_join(transcripts, by = c("transcript_id" = "name"), multiple = "all") %>%
      rename(PTC = is_ptc) %>%
      select(transcript_id, ref_transcript_name, everything()) %>%
      mutate(
        log2fold = round(log2fold, 2),
        padj = padj %>% scientific()
      ) %>%
      select(-c(seqnames, start, end, width, strand, cdna_thick, cdna_blocks)) %>%
      distinct(transcript_id, contrasts, .keep_all = TRUE)

    output$table_transcript <- renderReactable({
      reactable(
        df %>% select(transcript_id, ref_transcript_name, cds_position, contrasts, everything()),
        defaultSorted = c("PTC"),
        defaultPageSize = 5,
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
          cds_position = colDef(
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
  const cds_source = cellInfo.row['source'];
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
  const cc = cellInfo.row['match'];
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
            width = 220,
            show = TRUE,
            html = TRUE,
            vAlign = "center",
            cell = JS("
function(cellInfo) {
  const kd = cellInfo.row['Knockdown'];
  const cl = cellInfo.row['cellline'];
  const ko = cellInfo.row['Knockout'] || 'NoKO'
  const clone = cellInfo.row['clone'] || 'NA'


  if (kd === null) {
    return 'Not tested';
  } else {
    return `<div><strong>${kd}</strong> <br><i>Cell line</i>: ${cl}
<br><i>Clone</i>: ${clone}; <i> KO</i>: ${ko}</small></div>`;
  }
}")
          ),
          padj = colDef(
            header = with_tooltip(
              "padj", "DTE adjusted p-value."
            ),
            filterable = FALSE,
            vAlign = "center",
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
            vAlign = "center",
            format = colFormat(digits = 2),
            filterable = FALSE
          ),
          Knockdown = colDef(
            show = FALSE
          ),
          source = colDef(
            show = FALSE
          ),
          cds_trackhub_url = colDef(
            show = FALSE
          ),
          cds_id = colDef(
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
          name = colDef(
            show = FALSE
          ),
          clone = colDef(
            show = FALSE
          ),
          thick = colDef(
            show = FALSE
          ),
          match = colDef(
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
          ptc = colDef(
            show = FALSE
          ),
          itemRgb = colDef(
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
        p1 <- plot_annotation_cdna(transcripts)
        y_labs <- lapply(
          ggplot_build(p1)$layout$panel_params,
          \(x) x$y$get_labels()
        ) |>
          unlist()
        validate(need(length(y_labs) < 10,
                      "Too many transcripts, download plot instead."))
        p2 <- build_dotplot(df, y_labs)
        p_final <- p1 + p2 +
          patchwork::plot_layout(widths = c(4, 1), guides = "collec") &
          theme(
            legend.box = "vertical",
            legend.position = "left",
            legend.margin = margin(),
            legend.text = element_text(size = 6),
            legend.title = ggplot2::element_text(size = 6, face = "bold"),
            legend.key.size = unit(0.9, "line"),
            legend.key.width = unit(1, "line"),
            legend.key.height = unit(0.4, "line"),
            panel.spacing.x = unit(0.0, "cm")
          )
        ggsave(fname)
        p_final
      },
      res = 150
    )

    output$downloadPlot <- downloadHandler(
      filename = function() {
        basename(fname)
      },
      content = function(file) {
        file.copy(fname, file)
      }, contentType = "image/png"
    )
  })
}
