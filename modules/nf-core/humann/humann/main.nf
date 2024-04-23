process HUMANN_HUMANN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0':
        'biocontainers/humann:3.8--pyh7cba7a3_0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(metaphlan_profile)
    path chocophlan_db
    path uniref_db

    output:
    tuple val(meta), path("*_genefamilies.tsv") , emit: genefamilies
    tuple val(meta), path("*_pathabundance.tsv"), emit: pathabundance
    tuple val(meta), path("*_pathcoverage.tsv") , emit: pathcoverage
    tuple val(meta), path("*.log")              , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    humann \\
        --input ${reads} \\
        --output ./ \\
        --threads ${task.cpus} \\
        --taxonomic-profile ${metaphlan_profile} \\
        --nucleotide-database ${chocophlan_db} \\
        --protein-database ${uniref_db} \\
        --o-log ${prefix}.log \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genefamilies.tsv
    touch ${prefix}_pathabundance.tsv
    touch ${prefix}_pathcoverage.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
