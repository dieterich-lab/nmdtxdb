#
#
#
golem::fill_desc(
  pkg_name = "nmdtx",
  pkg_title = "NMD transcriptome",
  pkg_description = "PKG_DESC.",
  author_first_name = "Thiago",
  author_last_name = "Britto-Borges",
  author_email = "thiago.brittoborges@uni-heildeberg.de",
  repo_url = NULL
)
golem::set_golem_options()
usethis::use_mit_license("Thiago Britto-Borges")
usethis::use_readme_rmd(open = FALSE)
usethis::use_code_of_conduct(contact = "Thiago Britto-Borges")
usethis::use_lifecycle_badge("Experimental")
usethis::use_news_md(open = FALSE)
usethis::use_git()
golem::use_recommended_tests()
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)
rstudioapi::navigateToFile("dev/02_dev.R")
