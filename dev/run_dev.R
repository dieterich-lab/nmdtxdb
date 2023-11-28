options(golem.app.prod = FALSE, warn = 1)
options(shiny.fullstacktrace = TRUE)
golem::detach_all_attached()
golem::document_and_reload()
run_app(launch.browser = TRUE)

