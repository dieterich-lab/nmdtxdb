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
mod_transcript_ui <- function(id) {
  ns <- NS(id)

  reactableOutput(ns("table_transcript")) %>%
    shinycssloaders::withSpinner()
}


transcript_plot <- function(gtf) {
  if (nrow(gtf) < 1){ return("No CDS source to show.") }

  exons <- gtf %>% dplyr::filter(type == "exon")
  cds <- gtf %>% dplyr::filter(type == "CDS")
  introns <- to_intron(exons, group_var = "Name")
  feat_colors <- c("TRUE" = "firebrick", "FALSE" = "black")

  p1 <- exons %>%
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
    ) +
    scale_fill_manual(values = feat_colors) +
    theme_minimal()
    labs(y = "") +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "top"
    )

  htmltools::plotTag(p1, alt = "plots")

}

#' transcript view Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom purrr map
#' @importFrom reactable getReactableState renderReactable
#' @import stringr
#' @noRd
mod_transcript_server <- function(id, conn, tx, contrast, cds) {
  moduleServer(id, function(input, output, session) {
    output$table_transcript <- renderReactable({

      dte <- tbl(conn, "dte") %>%
        filter(
          transcript_id %in% !!tx,
          contrasts %in% !!contrast
        ) %>%
        select(transcript_id, contrasts, padj, log2fold) %>%
        collect() %>%
        left_join(load_metadata(conn), by = "contrasts") %>%
        select(-name)

      l2fc <- tbl(conn, "dte") %>%
        filter(contrasts %in% !!contrast) %>%
        collect()

      anno <- tbl(conn, "anno") %>%
        filter(transcript_id %in% !!tx) %>%
        select(ref_gene_name, transcript_id, ref_transcript_name,
               ref_transcript_id, lr_support, class_code) %>%
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
            \(x)  create_trackhub_url(position = x)) %>% unlist()
        )

      not_transcrips <- tbl(conn, "gtf") %>%
        filter(transcript_id %in% !!tx, type != "transcript") %>%
        collect() %>%
        left_join(
          x = transcripts %>% select(Name, cds_source, PTC),
          y = select(., !c(cds_source)),
          by = "Name", multiple = "all")


      df <- anno %>%
        left_join(dte, by = "transcript_id", multiple = "all") %>%
        distinct() %>%
        left_join(transcripts, by = "transcript_id", multiple = "all") %>%
        select(transcript_id, ref_transcript_name, PTC, lr_support, everything()) %>%
        distinct() %>%
        mutate(
          log2fold = round(log2fold, 2),
          padj = padj %>% scales::scientific(),
          class_code = case_when(
            class_code == "=" ~ "Complete",
            class_code == "n" ~ "Retained introns",
            class_code %in% c("c", "j", "k") ~ "Splicing variants",
            TRUE ~ "Other"
          )
        )

      reactable(
        df %>% select(
          -c(seqnames, start, Name, end, color)),
        highlight = TRUE,
        wrap = FALSE,
        details = function(index) transcript_plot(
          not_transcrips %>% filter(
            transcript_id == df[[index, 'transcript_id']],
            cds_source == df[[index, 'cds_source']])),
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
          PTC = colDef(
            header = with_tooltip(
              "PTC", "\u2713 if the transcript has a PTC else \u274c."
            ),
            show = TRUE,
            width = 60,
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
              "transcript_id", "Internal assembly identifier."
            ),
            width = 160,
            show = TRUE,
            vAlign = "center",
          ),
          trackhub_url = colDef(
            header = with_tooltip(
              "trackhub_url", "Link to the UCSC Genome Browser Trackhub."
            ),
            width = 60,
            show = TRUE,
            align = "center",
            vAlign = "center",
            html = TRUE,
            cell = JS("
function(cellInfo) {
  const url = cellInfo.row['trackhub_url'];

  if (url === undefined) {
    return '';
  } else {
    return `<a href='${url}' target='_blank'>URL</a>`;
  }
}")
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

  if (ti === undefined) {
    return lab;
  } else {
    const url = 'http://www.ensembl.org/id/' + ti;
    return `<a href='${url}' target='_blank'>${lab1}</a>${lab2}`;
  }
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
  })
}
