
#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' see https://github.com/ColinFay/golemexamples/blob/master/golemfuture/app.R
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(onStart = NULL,
                    enableBookmarking = NULL,
                    uiPattern = "/",
                    ...) {
  options <- list(
    warn = 0,
    shiny.port = 3838,
    shiny.host = "0.0.0.0"
  )

  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
