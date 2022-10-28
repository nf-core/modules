process FQ_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::fq=0.9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.9.1--h9ee0642_0':
        'quay.io/biocontainers/fq:0.9.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)
    path names

    output:
    tuple val(meta), path(""), emit: fastq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    fq filter \\
        --names $names \\
        $fastq > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq filter --version | sed 's/fq-filter //g'))
    END_VERSIONS
    """
}
