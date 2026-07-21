process LAST_POSTMASK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/05ff1dba1dcf86a28dd9d911ffcbc925175b370e77b5ea699ae421b8481556bc/data'
        : 'community.wave.seqera.io/library/last_gzip:40e65bc1d7d8624e'}"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    set -o pipefail
    last-postmask $args $maf | gzip --no-name > ${prefix}.maf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo "" | gzip > ${prefix}.maf.gz
    """
}
