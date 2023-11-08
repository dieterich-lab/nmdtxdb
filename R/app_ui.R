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
          img(
            src = "www/logo.svg",
            style = "width: 100%;"),
          br(),
          uiOutput("gene_info"),
          br(),
          selectizeInput(
            inputId = "gene_select",
            label = h3("Pick a gene:"),
            choices = NULL,
            selected = NULL,
            multiple = FALSE,
            options = list(create = FALSE)
          ),
          selectizeInput(
            inputId = "contrast_select",
            label = h3("Contrasts:"),
            choices = NULL,
            multiple = TRUE,
            options = list(
              create = FALSE,
              plugins = list("remove_button")
            )
          ),
          selectizeInput(
            inputId = "cds_source_select",
            label = h3("CDS sources:"),
            choices = NULL,
            multiple = TRUE,
            options = list(
              create = FALSE,
              plugins = list("remove_button")
            )
          ),
          div(
            action_button(
              input_id = "action_button_feedback",
              onclick = "window.open(
              'https://forms.gle/ZnaCwzNpFDPUHeh27',
              '_blank')",
              icon = icon("comments"),
              label = "Feedback",
              style = "width: 100%; padding: 0; height: 40px"
            )
          ),
          br(),
          div(
            action_button(
              input_id = "action_button_tutorial",
              onclick = "alert('Not implemented');",
              icon = icon("question"),
              label = "Tutorial",
              style = "width: 100%; padding: 0; height: 40px"
            )
          ),
        ),
        main_panel(
          width = 10,
          tabset(
            tabs = list(
              list(
                menu = "Introduction",
                content = mod_intro_ui("intro_1"),
                id = "intro_tab"
              ),
              list(
                menu = "Gene expression",
                content = mod_gene_ui("mod_gene1"),
                id = "gene_view_tab"
              ),
              list(
                menu = "Transcript expression",
                content = mod_transcript_ui("mod_transcript1"),
                id = "transcript_view_tab"
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
    ),
  )
}
