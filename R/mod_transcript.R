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
mod_transcript_server <- function(id, conn, select) {
  moduleServer(id, function(input, output, session) {
    anno <- reactive({
      validate(need(select, message = "Waiting selection"))
      tbl(conn, "anno") %>%
        dplyr::filter(gene_name == !!select)
    })

    output$table_transcript <- renderReactable({
      anno() %>%
        select(transcript_id) %>%
        distinct() %>%
        left_join(tbl(conn, "dtu2"), by = "transcript_id") %>%
        select(transcript_id, contrasts, padj, log2fold) %>%
        filter(!is.na(transcript_id)) %>%
        left_join(
          tbl(conn, "gtf") %>% select("transcript_name", "transcript_id", "transcript_biotype"),
          by = c("transcript_id")
        ) %>%
        distinct() %>%
        collect() %>%
        mutate(contrasts = str_sub(contrasts, 5)) %>%
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
