process PAIRTOOLS_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pairtools=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5' :
        'quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5' }"

    input:
    tuple val(meta), path(bam)
    path chromsizes

    output:
    tuple val(meta), path("*.pairsam.gz")  , emit: pairsam
    tuple val(meta), path("*.pairsam.stat"), emit: stat
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools \\
        parse \\
        -c $chromsizes \\
        $args \\
        --output-stats ${prefix}.pairsam.stat \\
        -o ${prefix}.pairsam.gz \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
