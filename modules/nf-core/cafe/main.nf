process CAFE {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1abaff43e5320ef242c17d05367cb50bf1e256c0e4f9211e522e0c187b7ab4ba/data'
        : 'community.wave.seqera.io/library/cafe:5.1.0--4ba563d3b90c6a6e'}"

    input:
    tuple val(meta), path(infile)
    path(tree)

    output:
    tuple val(meta), path("${prefix}") , emit: cafe
    path("$prefix/*_count.tab") , emit: cafe_base_count
    path("$prefix/*.tre") , emit: cafe_significant_trees
    path("$prefix/*_report.cafe") , emit: cafe_report
    path("$prefix/*results.txt") , emit: cafe_results
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('cafe'), val('5.1.0'), emit: versions_cafe, topic: versions
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
