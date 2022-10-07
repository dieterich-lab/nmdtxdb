transcript_view_grid <- grid_template(
  default = list(
    areas = rbind(c("left", "right")),
    cols_width = c("30%", "auto"),
    rows_height = c("auto")
  )
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
    # left = reactableOutput(ns("transcript_table")) %>%
    #   shinycssloaders::withSpinner(),
    # right = plotOutput(ns("transcript_structure")) %>%
    #   shinycssloaders::withSpinner()
  )
}

#' transcript view Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom reactable reactableOutput getReactableState renderReactable
#' @import stringr
#' @noRd
mod_transcript_server <- function(id, conn, select) {
  moduleServer(id, function(input, output, session) {

  })
}
