# Updating the docker container and making a new module release

witty.er is a commercial tool from Illumina. The container provided for the witty.er nf-core module is not provided nor supported by Illumina. Updating the witty.er versions in the container and pushing the update to Dockerhub needs to be done manually.
NOTE: an updated version of Dockerfile than the official one in github is used to build nf-core/wittyer:0.5.2, which is inserted here ./Dockerfile.

1. Navigate to the witty.er github repository. - [witty.er](https://github.com/Illumina/witty.er)
2. Download the latest release.
   ```bash
   wget https://github.com/Illumina/witty.er/archive/refs/tags/<VERSION>.tar.gz
   ```
3. Uncompress the released package.
   ```bash
   tar -xvf <VERSION>.tar.gz
   ```
4. Change to the uncompressed directory.
5. Build docker image using provided Dockerfile.

   ```bash
   docker build -t wittyer:<VERSION> --platform linux/amd64 .
   ```

6. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

   ```bash
   docker tag wittyer:<VERSION> quay.io/nf-core/wittyer:<VERSION>
   docker push quay.io/nf-core/wittyer:<VERSION>
   ```

7. Make the image public.
