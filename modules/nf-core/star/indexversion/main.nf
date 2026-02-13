process STAR_INDEXVERSION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data'
        : 'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4'}"

    output:
    path ("*.txt"), emit: index_version
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed "s/STAR_//"'), emit: versions_star, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "index_version"
    """
    STAR --help | grep 'versionGenome' | sed -e 's/versionGenome[ \t]*//' | sed -e 's/ .*//'  > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "index_version"
    """
    STAR --help | grep 'versionGenome' | sed -e 's/versionGenome[ \t]*//' | sed -e 's/ .*//'  > ${prefix}.txt
    """
}
