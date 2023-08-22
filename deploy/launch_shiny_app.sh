#!/bin/bash

R -e "renv::install('bioc::GenomicRanges')"
R -e "renv::restore()"
R -e "remotes::install_local('/app.tar.gz',upgrade='never')"
rm /app.tar.gz

setsid R -f rstart.R &> /criu_dumps/rlogs.log 2>&1 &
while ! grep -q -F 'Library loaded' /criu_dumps/rlogs.log
do
  echo 'Waiting for library load ...'
  sleep 5
done
rpid=$(pidof R)
criu-ns dump -t $rpid -vvv -o /criu_dumps/dump.log -D /criu_dumps/ && echo OK
touch /criu_dumps/dump.done
read -p 'Waiting for input'