library(nmdtx)
options("golem.app.prod" = TRUE)
db <- readRDS("/data/database.RDS")
gc()

while (! file.exists('/restore.marker')) {
    print ('waiting....')
    print('Library loaded')
    Sys.sleep(0.5)
}
system(paste('mkdir -p', tempdir()))
system('tail --follow=name /app/rlogs.log > /app/rlogs_container.log &')
shiny_user  <- basename(readLines('/shiny_user.txt'));
shiny_groups  <- basename(readLines('/shiny_usergroups.txt'));
print(paste('User', shiny_user))
Sys.setenv('SHINYPROXY_USERNAME' = shiny_user, 'SHINYPROXY_USERGROUPS' = shiny_groups)
print(Sys.getenv('SHINYPROXY_USERNAME'))
print(Sys.getenv('SHINYPROXY_USERGROUPS'))
options('shiny.port'=3838, shiny.host='0.0.0.0')
nmdtx::run_app()