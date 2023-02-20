# Updating the docker container and making a new module release

bcl2fastq is a commercial tool from Illumina. The container provided for the bcl2fastq nf-core module is not provided nor supported by Illumina. Updating the bcl2fastq versions in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download page. - [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/downloads.html): download the rpm of the desired bcl2fastq version with `curl` or `wget`.
2. Unpack the RPM package using `rpm2cpio bcl2fastq2-*.rpm | cpio -i --make-directories`.
   - Move the executable located in `<unpack_dir>/usr/bin/bcl2fastq` into the same folder as the Dockerfile.
   - Move the `css` and `xsl` directories from `<unpack_dir>/local/share/` into the same folder as the Dockerfile
3. Create and test the container:

   ```bash
   docker build . -t nfcore/bcl2fastq:<VERSION>
   ```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

   ```bash
   docker push nfcore/bcl2fastq:<VERSION>
   ```
