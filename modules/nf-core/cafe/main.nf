process CAFE {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cafe:5.1.0--h43eeafb_0':
        'biocontainers/cafe:5.1.0--h43eeafb_0' }"

    input:
    path(infile)
    path(tree)

    output:
    path(".")     , emit: cafe
    path("versions.yml")   , emit: versions
    path("*_count.tab") , emit: cafe_base_count
    path("*.tre") , emit: cafe_significant_trees
    path("*_report.cafe") , emit: cafereport
    path("*results.txt") , emit: cafe_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tr '\\r' '\\n' < $infile > infile.txt
    tr '\\r' '\\n' < $tree > treefile.txt
    cafe5 \\
        -i infile.txt \\
        -t treefile.txt \\
        $args \\
        --cores ${task.cpus} \\
        -o .


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch .
    touch versions.yml
    touch *_count.tab
    touch *.tre
    touch *_report.cafe
    touch *results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """
}
