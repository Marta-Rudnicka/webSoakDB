#!/bin/bash
# This script uses podman instead of Docker because that's the RedHat version.
module load pollux # Diamonds Kubernetes deployment (will require logging in once a week)
# podman image prune -a
podman build --rm -t web_soak_db .
x=$(podman images | grep -e 'web_soak_db' | awk '{print $3}')
podman tag $x gcr.io/diamond-privreg/xchemapps/web_soak_db:latest
podman push gcr.io/diamond-privreg/xchemapps/web_soak_db:latest
echo 'Finished building Container'
