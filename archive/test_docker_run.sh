#!/usr/bin/env bash

cd /home/romunov
echo "i'm here" > output.txt

echo "changing folders"
cd /home/romunov/biotools

echo "running docker"
docker run -dt --name name_legacy -v /home/romunov/:/mydata resolwebio/legacy:1.0.0 bash --login
echo "executing command"
docker exec -w /mydata name_legacy ls | cat > output.txt

echo "shutting down docker"
docker stop name_legacy
docker rm name_legacy