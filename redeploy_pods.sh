#!/bin/bash
module load pollux
x=$(kubectl get pods | grep -e 'webSoakDB' | awk '{print $1}')
kubectl delete pod $x
# Print new pods
kubectl get pods