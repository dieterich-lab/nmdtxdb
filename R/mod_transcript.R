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
mod_transcript_server <- function(id, conn, t_name, contrast) {
  moduleServer(id, function(input, output, session) {
    output$table_transcript <- renderReactable({
      conn %>%
        tbl("dtu2") %>%
        filter(
          transcript_name %in% !!t_name,
          contrasts %in% !!contrast
        ) %>%
        select(transcript_name, contrasts, transcript_biotype, padj, log2fold) %>%
        distinct() %>%
        collect() %>%
        reactable(
          .,
          language = reactableLang(
            filterPlaceholder = "Filter"
          ),
          filterable = TRUE,
          striped = TRUE,
          defaultSorted = c("padj"),
          showPageSizeOptions = TRUE,
          defaultPageSize = 5,
          pageSizeOptions = c(5, 10, 25, 50),
          highlight = TRUE,
          wrap = FALSE,
          rowStyle = list(cursor = "pointer"),
          theme = reactableTheme(
            stripedColor = "#f6f8fa",
            highlightColor = "#f0f5f9",
            cellPadding = "8px 12px",
            rowSelectedStyle = list(
              backgroundColor = "#eee",
              boxShadow = "inset 2px 0 0 0 #FF0000"
            )
          ),
          defaultColDef = colDef(
            sortNALast = TRUE
          ),
          columns = list(
            transcript_id = colDef(
              show = FALSE
            ),
            padj = colDef(
              format = colFormat(digits = 2),
              filterable = FALSE
            ),
            log2fold = colDef(
              name = "log2fc",
              format = colFormat(digits = 2),
              filterable = FALSE
            ),
            contrasts = colDef(
              width = 200
            )
          )
        )
    })
  })
}
