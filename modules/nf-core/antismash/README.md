# Updating the docker container and making a new module release

antiSMASH is an open source tool.
The developers do provide two premade docker containers, one that contains a very large database and one that does not.
The database containing version takes an extremely long time to pull as it's multiple gigabytes in size.
While theoretically we could use the one without the database, and we would prefer to use a container that is sync with the conda recipe.

Unfortunately, the default biocontainer recipe does not work out of the box, as the tool requires caching information to be stored inside the installation directories of the tool.
This caching procedure happens after database downloading, however the database does not fit on the bioconda/biocontainer CI nodes, thus the caching step cannot be executed.
Therefore, when supplying an externally downloaded database to the container, the tool tries to update the cache locations inside the installation directories which thus fails as containers are read only.

Therefore, for docker and singularity containers, we need to create a new container based on the biocontainers container, however we locally download the database, generate the cache within the container, _but_ then remove the database before finishing the build.

_Thanks to @mberacochea for finding the solution and providing the Dockerfile below!_

Updating the antiSMASH version in the container and pushing the update to Dockerhub needs to be done manually.

1. Ensure the user that will push the container has access to the nf-core organization on Dockerhub.
2. Change into the `modules/nf-core/antismash` directory:

   ```bash
   cd modules/nf-core/antismash
   ```

3. Create and test the container:

   ```bash
   docker build . -t quay.io/nf-core/antismash:<version>_<build>
   ```

4. Push the container to Quay:

   ```bash
   docker push quay.io/nf-core/antismash:<version>_<build>
   ```

5. On [https://quay.io/organization/nf-core](https://quay.io/organization/nf-core), select the new repository, press actions and make public.
