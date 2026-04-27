process BAMUTIL_SQUEEZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamutil:1.0.15--h5ca1c30_6':
        'biocontainers/bamutil:1.0.15--h5ca1c30_6' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('bamutil'), eval("bam help 2>&1 | sed -n 's/.*Version: \([^;]*\).*/\1/p'"), emit: versions_bamutil, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_squeezed"
    def ref = fasta ? "--refFile ${fasta}" : ''
    """
    bam squeeze \
        --in $bam \
        --out ${prefix}.bam \
        ${ref} \
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_squeezed"
    """
    touch ${prefix}.bam
    """
}
