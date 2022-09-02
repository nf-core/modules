process HMMER_ESLREFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1':
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    tuple val(meta), val(format), path(alignment)

    output:
    tuple val(meta), path("*.aln.gz"), emit: aln
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    esl-reformat \\
        $format \\
        $alignment \\
        $args \\
        -o ${prefix}.aln \\
        -T $prefix \\
        $bam

    gzip ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
