#' intro UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import markdown
mod_intro_ui <- function(id){
  ns <- NS(id)
  this_grid <- create_grid(
    areas = rbind(c("left")),
    cols_width = c("auto"),
    rows_height = c("auto")
  )

  grid(
    this_grid,
    left = includeMarkdown(app_sys("app/www/intro.md"))
  )
}

#' intro Server Functions
#'
#' @noRd
mod_intro_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}
