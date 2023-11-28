aout <- pkgbuild::build(
  path = ".",
  dest_path = 'deploy',
  vignettes = FALSE,
)

system("scp database2.RDS tbrittoborges@shinynew:apps/nmd-app/data/database.RDS")
system("scp deploy/nmdtx_0.0.0.9000.tar.gz tbrittoborges@shinynew:apps/nmd-app/deploy/")
# ssh shinynew
# cd apps/nmd-app/deploy/
# ./deploy/build.sh
# sudo -i
# cd /tmp_shiny/nmd_transcriptome_criu_dump/
# rm *
# exit
# sudo mv criu_dumps/* /tmp_shiny/nmd_transcriptome_criu_dump/
# sudo rmate /srv/shiny/prod/application.yml
# sudo systemctl restart shinyproxy-prod
#

