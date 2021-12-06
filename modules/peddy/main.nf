process PEDDY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::peddy=0.4.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/peddy:0.4.8--pyh5e36f6f_0' :
        'quay.io/biocontainers/peddy:0.4.8--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path ped

    output:
    tuple val(meta), path("*.html")     , emit: html
    tuple val(meta), path("*.csv")      , emit: csv
    tuple val(meta), path("*.peddy.ped"), emit: ped
    tuple val(meta), path("*.png")      , emit: png
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    peddy \\
        $args \\
        --plot \\
        -p $task.cpus \\
        $vcf \\
        $ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peddy: \$( peddy --version 2>&1 | sed 's/peddy, version //' )
    END_VERSIONS
    """
}
