process PRESEQ_CCURVE {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::preseq=3.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.1.2--h445547b_2':
        'quay.io/biocontainers/preseq:3.1.2--h445547b_2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.c_curve.txt"), emit: c_curve
    tuple val(meta), path("*.log")        , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        c_curve \\
        $args \\
        $paired_end \\
        -output ${prefix}.c_curve.txt \\
        $bam
    cp .command.err ${prefix}.command.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
}
