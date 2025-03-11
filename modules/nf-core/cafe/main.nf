process CAFE {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cafe:5.1.0--h43eeafb_0':
        'biocontainers/cafe:5.1.0--h43eeafb_0' }"

    input:
    tuple val(meta), path(infile)
    path(tree)

    output:
    tuple val(meta), path("${prefix}") , emit: cafe
    path("versions.yml") , emit: versions
    path("$prefix/*_count.tab") , emit: cafe_base_count
    path("$prefix/*.tre") , emit: cafe_significant_trees
    path("$prefix/*_report.cafe") , emit: cafe_report
    path("$prefix/*results.txt") , emit: cafe_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tr '\\r' '\\n' < $infile > infile.txt
    tr '\\r' '\\n' < $tree > treefile.txt
    cafe5 \\
        -i infile.txt \\
        -t treefile.txt \\
        $args \\
        --cores ${task.cpus} \\
        -o ${prefix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir ${prefix}
    touch versions.yml
    touch ${prefix}/*_count.tab
    touch ${prefix}/*.tre
    touch ${prefix}/*_report.cafe
    touch ${prefix}/*results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """
}
