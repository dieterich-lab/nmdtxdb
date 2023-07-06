transcript_view_grid <- create_grid(
  rbind(c("top")),
  c("auto"),
  c("auto")
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
mod_transcript_ui <- function(id) {
  ns <- NS(id)

  grid(
    transcript_view_grid,
    top = reactableOutput(ns("table_transcript")) %>% shinycssloaders::withSpinner()
  )
}

#' transcript view Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom reactable getReactableState renderReactable
#' @import stringr
#' @noRd
mod_transcript_server <- function(id, conn, tx, contrast, cds) {
  moduleServer(id, function(input, output, session) {
    output$table_transcript <- renderReactable({
      dte <- conn %>%
        tbl("dte") %>%
        filter(
          transcript_id %in% !!tx,
          contrasts %in% !!contrast
        ) %>%
        select(transcript_id, contrasts, padj, log2fold)


      gtf <- conn %>%
        tbl("gtf") %>%
        filter(
          transcript_id %in% !!tx,
          type == "transcript",
          cds_source %in% !!cds
        ) %>%
        select(transcript_id, cds_source, lr_support, color, class_code) %>%
        mutate(PTC = as.character(color == "#FF0000")) %>%
        select(-color)

      gtf %>%
        left_join(dte, by = "transcript_id") %>%
        collect() %>%
        mutate(
          log2fold = round(log2fold, 2),
          padj = padj %>% scales::scientific()
        ) %>%
        reactable(
          .,
          language = reactableLang(
            filterPlaceholder = "Filter"
          ),
          filterable = TRUE,
          striped = TRUE,
          defaultSorted = c("padj"),
          showPageSizeOptions = TRUE,
          defaultPageSize = 10,
          pageSizeOptions = c(5, 10, 25, 50),
          highlight = TRUE,
          wrap = FALSE,
          theme = reactableTheme(
            highlightColor = "#FFFFBF",
            cellPadding = "8px 12px"
          ),
          defaultColDef = colDef(
            sortNALast = TRUE
          ),
          columns =
            list(
              transcript_id = colDef(
                width = 180
              ),
              contrasts = colDef(
                width = 200
              )
            )
          #   padj = colDef(
          #     format = colFormat(digits = 2),
          #     filterable = FALSE
          #   ),
          #   log2fold = colDef(
          #     name = "log2fc",
          #     format = colFormat(digits = 2),
          #     filterable = FALSE
          #   ),
          #   contrasts = colDef(
          #     width = 200
          #   )
          # )
        )
    })
  })
}
