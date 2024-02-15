process CENTRIFUGE_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4.1--hdcf5f25_1' :
        'biocontainers/centrifuge:1.0.4.1--hdcf5f25_1' }"

    input:
    tuple val(meta), path(fasta)
    path conversion_table
    path taxonomy_tree
    path name_table
    path size_table

    output:
    tuple val(meta), path("*.cf") , emit: cf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def size_table_cmd = size_table ? "--size_table ${size_table}" : ""
    """
    centrifuge-build \\
        -p $task.cpus \\
        $fasta \\
        ${prefix} \\
        --conversion-table $conversion_table \\
        --taxonomy-tree $taxonomy_tree \\
        --name-table $name_table \\
        ${size_table_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.cf
    touch ${prefix}.2.cf
    touch ${prefix}.3.cf
    touch ${prefix}.4.cf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}
