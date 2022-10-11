#
#
#
devtools::check()
rhub::check_for_cran()
devtools::build()
golem::add_dockerfile_with_renv_shinyproxy()
