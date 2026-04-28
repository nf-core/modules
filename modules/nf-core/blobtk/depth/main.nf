process BLOBTK_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
<<<<<<< HEAD
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a04d3e2b8c1b50d6e26463114970a68afeddb5a0f302177cf1ed6b8b1f82b1c/data' :
        'community.wave.seqera.io/library/pip_blobtk:2fbb8045b9880daf' }"
=======
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/243d043f1c9e152e75dbb0ef8c64022df50efbcaa4e1bbaea36bebd751e84e93/data' :
        'community.wave.seqera.io/library/blobtk:0.7.1--e3f63bb2cdc8fb96' }"
>>>>>>> master

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
