process FCSGX_CLEANGENOME {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    tuple val(meta), path(fasta), path(fcsgx_report)

    output:
    tuple val(meta), path("*.cleaned.fasta")     , emit: cleaned
    tuple val(meta), path("*.contaminants.fasta"), emit: contaminants
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gx \\
        clean-genome \\
        --input ${fasta} \\
        --action-report ${fcsgx_report} \\
        --output ${prefix}.cleaned.fasta \\
        --contam-fasta-out ${prefix}.contaminants.fasta \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cleaned.fasta
    touch ${prefix}.contaminants.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """
}
