process GAPPA_EXAMINEHEATTREE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gappa:0.8.0--h9a82719_0':
        'biocontainers/gappa:0.8.0--h9a82719_0' }"

    input:
    tuple val(meta), path(jplace)

    output:
    tuple val(meta), path("*.newick")     , emit: newick  , optional: true
    tuple val(meta), path("*.nexus")      , emit: nexus   , optional: true
    tuple val(meta), path("*.phyloxml")   , emit: phyloxml, optional: true
    tuple val(meta), path("*.svg")        , emit: svg     , optional: true
    tuple val(meta), path("*.colours.txt"), emit: colours
    tuple val(meta), path("*.log")        , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gappa \\
        examine \\
        heat-tree \\
        --threads $task.cpus \\
        --file-prefix ${prefix}. \\
        --jplace-path $jplace \\
        $args \\
        --log-file ${prefix}.log

    grep '^ *At' ${prefix}.log > ${prefix}.colours.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gappa: \$(echo \$(gappa --version 2>&1 | sed 's/v//' ))
    END_VERSIONS
    """
}
