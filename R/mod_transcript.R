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

      gtf <- tbl(conn, "gtf") %>%
        filter(
          transcript_id %in% !!tx,
          type == "transcript",
          cds_source %in% !!cds
        ) %>%
        select(transcript_id, cds_source, lr_support, color, class_code, seqnames, start, end) %>%
        mutate(PTC = as.character(color == "#FF0000")) %>%
        select(-color) %>%
        collect() %>%
        mutate(
          position = position_from_gtf(.),
          trackhub_url = purrr::map(position, \(x)  create_trackhub_url(position = x)) %>% unlist()
        )


      anno <- tbl(conn, "anno") %>%
        filter(transcript_id %in% !!tx) %>%
        select(transcript_id, ref_transcript_name, ref_transcript_id) %>%
        collect()

      gtf <- gtf %>%
        left_join(dte, by = "transcript_id") %>%
        distinct() %>%
        left_join(anno, by = "transcript_id") %>%
        collect() %>%
        select(transcript_id, ref_transcript_name, class_code, everything()) %>%
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
        gtf,
        highlight = TRUE,
        wrap = FALSE,
        theme = reactableTheme(
          borderColor = "#dfe2e5",
          stripedColor = "#f6f8fa",
          highlightColor = "#FFFFBF",
          cellPadding = "8px 12px"
        ),
        defaultColDef = colDef(
          sortNALast = TRUE,
          show = FALSE
        ),
        columns = list(
          PTC = colDef(
            name = "PTC",
            show = TRUE,
            width = 60,
            align = "center",
            vAlign = "center",
            cell = function(value) {
              if (value == "false") "\u274c" else "\u2713"
            }
          ),
          lr_support = colDef(
            name = "LRS",
            width = 60,
            show = TRUE,
            align = "center",
            vAlign = "center",
            cell = function(value) {
              if (value == "FALSE") "\u274c" else "\u2713"
            }
          ),
          transcript_id = colDef(
            name = "Transcript id",
            width = 160,
            show = TRUE,
            vAlign = "center",
          ),
          trackhub_url = colDef(
            name = "TrackHub",
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
            name = "Ref Transcript",
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
            width = 160,
            show = TRUE,
            html = TRUE,
            vAlign = "center",
            cell = JS("
    function(cellInfo) {
    const kd = cellInfo.row['Knockdown']
    const cl = cellInfo.row['cellline']
    const ko = cellInfo.row['Knockout'] || 'NoKO'
        return (
          '<div>' +
          '<strong>' + kd + '</strong> <br>' +
          '<small><i>Cell-line</i>: ' + cl +
          ';<i> KO</i>: ' + ko + ' </small></div>'
        )
      }")
          ),
          padj = colDef(
            filterable = FALSE,
            show = TRUE,
            width = 100,
            align = "right"
          ),
          log2fold = colDef(
            name = "l2fc",
            show = TRUE,
            width = 80,
            format = colFormat(digits = 2),
            filterable = FALSE
          )
        )
      )
    })
  })
}
