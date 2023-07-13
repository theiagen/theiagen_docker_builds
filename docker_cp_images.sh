#!/bin/bash
# This script will take a file containing the paths to docker images, such as ncbi/sra-human-scrubber:2.1.0 or quay.io/ncbi/sra-human-scrubber:2.1.0
# then will pull them locally, and push them up to the theiagen Google Artifact Registry
# https://console.cloud.google.com/artifacts/browse/general-theiagen?orgonly=true&project=general-theiagen&supportedpurview=organizationId
# You need to be logged in on the command line with gcloud
# and for each directory (like ncbi) you need to use the web interfact to create it (and set read permissions to alUsers).

file_of_dockers=$1

images=$(cat file_of_dockers)
AR_PROJECT=general-theiagen

if [ -z "${AR_PROJECT}" ]
then
    echo ERROR: AR_PROJECT must be set before running this
    exit 1
fi

for img in ${images}
do
    dest_img=${img}
    # if 'quay.io/' is at the start of the img name, remove it and store in a new variable
    if [[ ${img} == quay.io/* ]]
    then
        dest_img=$(echo ${img} | sed 's/quay.io\///')
    else
        dest_img=${img}
    fi
    echo us-docker.pkg.dev/${AR_PROJECT}/${dest_img}
    docker pull ${img}
    docker tag ${img} us-docker.pkg.dev/${AR_PROJECT}/${dest_img}
    docker push us-docker.pkg.dev/${AR_PROJECT}/${dest_img}
done
