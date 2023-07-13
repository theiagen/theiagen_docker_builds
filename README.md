# Theiagen hosted docker images


## Production (public)
Use the docker repository for production images. This is a public repository, so anyone can pull from it. This is where we will push images that are used for production.
Images in this repository should not be deleted.

```bash
# download code
git clone https://github.com/theiagen/theiagen_docker_builds.git

# Login to Google Cloud
gcloud auth configure-docker us-docker.pkg.dev

# Build the docker image locally for GCP
docker build -t us-docker.pkg.dev/general-theiagen/theiagen/myapp:1.2.3 myapp:1.2.3

# Push the docker image to GCP
docker push us-docker.pkg.dev/general-theiagen/theiagen/myapp:1.2.3
```

To use the image locally:
```bash
docker pull us-docker.pkg.dev/general-theiagen/theiagen/myapp:1.2.3
```

## Private
Use the docker-private repository for private images. You need to be authenticated to pull from this repository. It is where you can push images which may contain information which needs to remain private, for example where private information must be included in a build or where data/software licences do not allow for redistribution. Usage of this repository should be kept to a minimum.

```bash
# download code
git clone https://github.com/theiagen/theiagen_docker_builds.git

# Login to Google Cloud
gcloud auth configure-docker us-docker.pkg.dev

# Build the docker image locally for GCP
docker build -t us-docker.pkg.dev/general-theiagen/docker-private/myapp:1.2.3 myapp:1.2.3

# Push the docker image to GCP
docker push us-docker.pkg.dev/general-theiagen/docker-private/myapp:1.2.3
```

To use the image:
```bash
gcloud auth configure-docker us-docker.pkg.dev
docker pull us-docker.pkg.dev/general-theiagen/docker-private/myapp:1.2.3
```
# Migrating docker images from other services
Use the docker_cp_images.sh script to pull images from other services like quay.io or dockerhub and push them to theiagens google artifact registry. If there is a new organisation (directory), you should create a new repository on the Google Artifact Registry first and set the permissions to be read for allUsers. Do not allow public writes.  Set the location to be multi-region (us) as of the images will be used by VMs running in that region, and network charges are zero.

```bash
# download code
git clone https://github.com/theiagen/theiagen_docker_builds.git

# create a file with 1 image name and tag per line (repo/myapp:1.2.3, or quay.io/repo/myapp:1.2.3)
bash docker_cp_images.sh <image_list.txt>
```


# Depricated - Public repositories (not hosted by Theiagen)
To build the docker image locally and push to dockerhub and quay, do the following. This will require write access to Theiagen's dockerhub and quay organizations.

Theiagen Quay repos: https://quay.io/organization/theiagen

Theiagen Dockerhub repos: https://hub.docker.com/u/theiagen

```bash
# download code
git clone https://github.com/theiagen/theiagen_docker_builds.git

# build docker image locally, for dockerhub
docker build -t theiagen/terra-tools:2023-07-05 broad-terra-tools-plus/2023-07-05/

# build docker image locally, but instead for quay.io
docker build -t quay.io/theiagen/terra-tools:2023-07-05 broad-terra-tools-plus/2023-07-05/

# push images to dockerhub and quay
docker push theiagen/terra-tools:2023-07-05
docker push quay.io/theiagen/terra-tools:2023-07-05

```