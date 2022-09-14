process HMMER_ESLREFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1':
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    tuple val(meta), path(seqfile)

    output:
    tuple val(meta), path("*.sequences.gz"), emit: seqreformated
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    esl-reformat \\
        -o ${prefix}.sequences \\
        $args \\
        $seqfile

    gzip ${prefix}.sequences

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer/easel: \$(esl-reformat -h | grep -o '^# Easel [0-9.]*' | sed 's/^# Easel *//')
    END_VERSIONS
    """
}
