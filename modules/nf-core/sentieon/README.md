# Updating the docker container and making a new module release

Sentieon is a commercial tool.
The container provided for the nf-core modules is not provided by Sentieon.
But follow the example provided by Sentieon here: https://github.com/Sentieon/sentieon-docker
Updating the Sentieon versions in the container and pushing the update to quay.io needs to be done manually.

1. Create and test the container:

   ```bash
   docker build . -t quay.io/nf-core/sentieon:<VERSION> --build-arg SENTIEON_VERSION=<VERSION>
   ```

2. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

   ```bash
   docker push quay.io/nf-core/sentieon:<VERSION>
   ```
