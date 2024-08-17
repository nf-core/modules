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
    path("Out_folder")     , emit: cafe
    path("versions.yml")   , emit: versions
    // I want to add the below out paths, but when I do this. nf-test fails for snapshot. We cannot snapshot the entire out folder, as it has time stamps in files and random numbers
    //path("Out_folder/Base_count.tab") , emit: cafe_base_count
    //path("Out_folder/*.tre") , emit: cafe_significant_trees
    //path("Out_folder/*_report.cafe") , emit: cafereport
    //path("Out_folder/*results.txt") , emit: cafe_results

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
        -o Out_folder


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch Out_folder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """
}
