aout <- pkgbuild::build(
  path = ".",
  dest_path = 'deploy',
  vignettes = FALSE,
)

system("scp database2.RDS tbrittoborges@shinynew:apps/nmd-app/data/database.RDS")
system("scp deploy/nmdtx_0.0.0.9000.tar.gz tbrittoborges@shinynew:apps/nmd-app/deploy/")
# ssh shinynew# cd apps/nmdtx/deploy/
# docker build -f Dockerfile --progress=plain -t nmdtx:dev .
#

