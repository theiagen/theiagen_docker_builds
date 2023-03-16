# theiagen_docker_builds

To build the docker image locally and push to dockerhub and quay, do the following. This will require write access to Theiagen's dockerhub and quay organizations.

```bash
# download code
git clone https://github.com/theiagen/theiagen_docker_builds.git

# build docker image locally, for dockerhub
docker build -t theiagen/terra-tools:2023-03-16 broad-terra-tools-plus/2023-03-16/

# build docker image locally, but instead for quay.io
docker build -t quay.io/theiagen/terra-tools:2023-03-16 broad-terra-tools-plus/2023-03-16/

# push images to dockerhub and quay
docker push theiagen/terra-tools:2023-03-16
docker push quay.io/theiagen/terra-tools:2023-03-16

```
