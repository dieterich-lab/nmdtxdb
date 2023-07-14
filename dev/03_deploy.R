devtools::check()
#rhub::check_for_cran()
devtools::build()
golem::add_dockerfile_with_renv_shinyproxy(output_dir = 'deploy/')
out <- pkgbuild::build(
  path = ".",
  dest_path = 'deploy',
  vignettes = FALSE,
)
# cd /repos/nmdtx/
# scp -r deploy/nmdtx_0.0.0.9000.tar.gz  tbrittoborges@shiny:apps/nmdtx/deploy/nmdtx_0.0.0.9000.tar.gz
#
# ssh shiny
# screen -S nmdtx:dev
# cd apps/nmdtx/deploy/
# docker build -f Dockerfile --progress=plain -t nmdtx:dev .
#

