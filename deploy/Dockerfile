FROM shinyverse_criu_base:latest

# setup app

COPY .Renviron .Renviron
COPY deploy/nmdtx_*.tar.gz /app.tar.gz
COPY deploy/rstart.R rstart.R
COPY deploy/restore.sh restore.sh
COPY deploy/launch_shiny_app.sh launch_shiny_app.sh
RUN mkdir -p /criu_dumps \
    && chmod 755 restore.sh launch_shiny_app.sh
EXPOSE 3838
CMD /bin/bash /restore.sh