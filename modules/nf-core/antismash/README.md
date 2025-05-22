# Updating the docker container and making a new module release

antiSMASH is an open source tool.
The developers do provide their own docker container but this uses some additional modifications that do not allow us to run out of the box.
Instead we prefer to use the biocontainers container for consistency with the bioconda recipe.
However, the tool requires a preparation caching step that depends on the location of installation files of the tool itself (which works with conda versions of the tool, as conda installation the database preparation can be executed again by the user and installation paths can be updated).

Therefore, for docker and singularity containers, we need to create a new contianer based on the biocontainers container, with the additional preparation step executed to generated the configuration _but_ then remove the database after.

_Thanks to @mberacochea for finding the solution and providing the Dockerfile below!_

Updating the antiSMASH version in the container and pushing the update to Dockerhub needs to be done manually.

1. Ensure the user that will push the container has access to the nf-core organization on Dockerhub.
2. Change into the `modules/nf-core/antismash` directory:

   ```bash
   cd modules/nf-core/antismash
   ```

3. Create and test the container:

   ```bash
   docker build . -t quay.io/nf-core/antismash:<TOOL>
   ```

4. Push the container to Dockerhub:

   ```bash
   docker push quay.io/nf-core/antismash:8.0.0
   ```

5. On [https://quay.io/organization/nf-core](https://quay.io/organization/nf-core), select the new repository, press actions and make public.
