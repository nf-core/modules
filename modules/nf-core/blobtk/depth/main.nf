process BLOBTK_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/08833d1b91f41024e06e2cb5a982598531199c04e6544885d42ef2cb0480de18/data' :
        'community.wave.seqera.io/library/blobtk:0.8.0--2fe0d833a26e0cd9' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path('*.regions.bed.gz') , emit: bed
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    blobtk depth \\
        -b ${bam} \\
        $args \\
        -O ${prefix}.regions.bed.gz \\
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.regions.bed.gz
    """
}
