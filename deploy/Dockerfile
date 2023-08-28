FROM shinyverse_criu_base:latest

WORKDIR /app

COPY .Renviron .Renviron
COPY renv.lock renv.lock
COPY deploy/nmdtx_*.tar.gz app.tar.gz
COPY deploy/rstart.R rstart.R
COPY deploy/restore.sh /restore.sh
COPY deploy/launch_shiny_app.sh launch_shiny_app.sh

RUN mkdir -p /criu_dumps \
    && chmod 755 /restore.sh launch_shiny_app.sh  \
    && mkdir -p /data 

EXPOSE 3838

CMD /bin/bash /restore.sh