aout <- pkgbuild::build(
  path = ".",
  dest_path = 'deploy',
  vignettes = FALSE,
)

install.packages(aout, repos = NULL, type ="source")
db <- readRDS("database.RDS")
devtools::test()

system("scp database.RDS tbrittoborges@shinynew:apps/nmd-app/data/database.RDS")
system(stringr::str_glue("scp {aout} tbrittoborges@shinynew:apps/nmd-app/deploy/"))
system("scp renv.lock tbrittoborges@shinynew:apps/nmd-app/")

# ssh shinynew
# cd apps/nmd-app/
# ./deploy/build.sh
# sudo -i
# cd /tmp_shiny/nmd_transcriptome_criu_dump/
# rm *
# exit
# sudo mv criu_dumps/* /tmp_shiny/nmd_transcriptome_criu_dump/
# sudo rmate /srv/shiny/prod/application.yml
# sudo systemctl restart shinyproxy-prod
#
# cp canonical.bb ensembl.bb riboseq.bb openprot.bb /prj/trackhubs/nmd_transcriptome/hg38/
# deploy_new_trackhub.sh /prj/trackhubs/nmd_transcriptome/
