process FGUMI_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e6/e613097ca7c84595a8683b6b042d113ac40e1936a809477aa95fcd1f8f3bfca2/data'
:         'community.wave.seqera.io/library/fgumi:0.6.0--c97194d17da0d1cd' }"

    input:
    tuple val(meta), path(bams, stageAs: "?/*")

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi \\
        merge \\
        --output ${prefix}.bam \\
        --threads $task.cpus \\
        ${args} \\
        ${bams}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
