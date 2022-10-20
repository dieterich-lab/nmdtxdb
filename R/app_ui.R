#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shiny.semantic
#' @importFrom reactable reactableOutput
#' @importFrom shinycssloaders withSpinner
#'
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    semanticPage(
      sidebar_layout(
        sidebar_panel(
          width = 2,
          selectizeInput(
            inputId = "gene_select",
            label = h2("Select a gene:"),
            choices = NULL,
            selected = NULL,
            multiple = FALSE,
            options = list(create = FALSE)
          ),
          uiOutput("gene_info"),
          div(a(
            href = "https://forms.gle/ZnaCwzNpFDPUHeh27",
            "Link to feedback form.",
            target = "_blank"
          )),
        ),
        main_panel(
          width = 10,
          tabset(
            tabs = list(
              list(
                menu = "Gene expression",
                content = mod_gene_ui("mod_gene1"),
                id = "gene_view_tab"
              ),
              list(
                menu = "Transcript table",
                content = mod_transcript_ui("mod_transcript1"),
                id = "transcript_view_tab"
              ),
              list(
                menu = "Advanced view",
                content = mod_phase1_ui("mod_phase1"),
                id = "advanced_tab"
              )
            ),
            active = "second_tab",
            id = "transcript_tabset"
          )
        )
      )
    )
  )
}
#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "nmdtx"
    )
  )
}
