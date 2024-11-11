process NANOQC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoqc:0.10.0--pyhdfd78af_0':
        'biocontainers/nanoqc:0.10.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*.log") , emit: log
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    nanoQC \\
        $fastq \\
        $args
    
    mv *.html ${prefix}.html
    mv *.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoqc: \$(nanoQC --version |& sed '1!d ; s/NanoQC //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}.html
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoqc: \$(nanoQC --version |& sed '1!d ; s/NanoQC //')
    END_VERSIONS
    """
}
