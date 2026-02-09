process PCGR_GETREF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/51/517cc3a46129fc586191412f29246889124f4654640a53230aa203a2bcc1d7dc/data'
        : 'community.wave.seqera.io/library/coreutils_curl_gzip_tar:17a6ea9a6766c02a'}"

    input:
    tuple val(meta), val(bundleversion), val(genome)

    output:
    tuple val(meta), path("${bundleversion}"), emit: pcgrref
    tuple val("${task.process}"), val('curl'),  eval("curl --version | head -1 | cut -d ' ' -f 2"), topic: versions, emit: versions_curl
    tuple val("${task.process}"), val('gzip'),  eval("gzip --version | head -1 | cut -d ' ' -f 2"), topic: versions, emit: versions_gzip
    tuple val("${task.process}"), val('tar'),   eval("tar --version | head -1 | cut -d ' ' -f 4"),  topic: versions, emit: versions_tar

    when:
    task.ext.when == null || task.ext.when

    script:
    // no args as this process is only for downloading and unpacking the reference data bundle, which doesn't require any additional arguments
    def bundle = "pcgr_ref_data.${bundleversion}.${genome}.tgz"
    """
    curl -O https://insilico.hpc.uio.no/pcgr/${bundle}
    gzip -dc ${bundle} | tar xvf -

    mkdir ${bundleversion}
    mv data/ ${bundleversion}
    """

    stub:
    """
    mkdir ${bundleversion}
    """
}
