process FREYJA_VARIANTS {
    tag "$meta.id"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.3--pyhdfd78af_0' :
        'biocontainers/freyja:2.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.variants.tsv"), path("*.depth.tsv"), emit: variants
    tuple val("${task.process}"), val('freyja'), eval("freyja --version | sed 's/.* //'"), topic: versions, emit: versions_freyja

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freyja \\
        variants \\
        $args \\
        --ref $fasta \\
        --variants ${prefix}.variants.tsv \\
        --depths ${prefix}.depth.tsv \\
        $bam

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.tsv
    touch ${prefix}.depth.tsv

    """
}
