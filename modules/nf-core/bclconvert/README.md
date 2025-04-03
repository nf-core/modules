# Updating the docker container and making a new module release

bcl-convert is a commercial tool from Illumina. The container provided for the bcl-convert nf-core module is not provided nor supported by Illumina. Updating the bcl-convert versions in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download page. - [BCL Convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert/downloads.html): find and download the latest version of bcl-convert RPM package.
2. Update the `BCL_VERSION` in the `Dockerfile`
3. Create and test the container:

   ```bash
   docker build . -t quay.io/nf-core/bclconvert:<VERSION>
   ```

4. Access rights are needed to push the container to the Dockerhub nf-core organization, please ask a core team member to do so.

   ```bash
   docker push quay.io/nf-core/bclconvert:<VERSION>
   ```
