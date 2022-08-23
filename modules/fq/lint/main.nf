process FQ_LINT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fq=0.9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/fq:0.9.1--h9ee0642_0':
        'quay.io/biocontainers/fq:0.9.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    path "*.fq_lint.txt", emit: lint
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fq lint \\
        $args \\
        $fastq > ${prefix}.fq_lint.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq lint --version | sed 's/fq-lint //g'))
    END_VERSIONS
    """
}
