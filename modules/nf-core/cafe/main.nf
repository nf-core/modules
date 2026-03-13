process CAFE {
    label 'process_high'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cafe:5.1.0--h43eeafb_0':
        'biocontainers/cafe:5.1.0--h43eeafb_0' }"

    input:
    tuple val(meta), path(infile)
    path(tree)

    output:
    tuple val(meta), path("${prefix}") , emit: cafe
    path("$prefix/*_count.tab") , emit: cafe_base_count
    path("$prefix/*.tre") , emit: cafe_significant_trees
    path("$prefix/*_report.cafe") , emit: cafe_report
    path("$prefix/*results.txt") , emit: cafe_results
    tuple val("${task.process}"), val('cafe'), eval('echo 5.1.0'), emit: versions_cafe, topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    tr '\\r' '\\n' < $infile > infile.txt
    tr '\\r' '\\n' < $tree > treefile.txt
    cafe5 \\
        -i infile.txt \\
        -t treefile.txt \\
        $args \\
        --cores ${task.cpus} \\
        -o ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir ${prefix}
    touch ${prefix}/*_count.tab
    touch ${prefix}/*.tre
    touch ${prefix}/*_report.cafe
    touch ${prefix}/*results.txt
    """
}
