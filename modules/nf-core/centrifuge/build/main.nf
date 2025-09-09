process CENTRIFUGE_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4.2--hdcf5f25_0'
        : 'biocontainers/centrifuge:1.0.4.2--hdcf5f25_0'}"

    input:
    tuple val(meta), path(fasta)
    path conversion_table
    path taxonomy_tree
    path name_table
    path size_table

    output:
    tuple val(meta), path("${prefix}/"), emit: cf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def size_table_cmd = size_table ? "--size_table ${size_table}" : ""
    """
    mkdir ${prefix}

    centrifuge-build \\
        -p ${task.cpus} \\
        ${fasta} \\
        ${prefix}/${prefix} \\
        --conversion-table ${conversion_table} \\
        --taxonomy-tree ${taxonomy_tree} \\
        --name-table ${name_table} \\
        ${size_table_cmd} \\
        ${args} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """

    stub:
    def _args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/
    touch ${prefix}/${prefix}.1.cf
    touch ${prefix}/${prefix}.2.cf
    touch ${prefix}/${prefix}.3.cf
    touch ${prefix}/${prefix}.4.cf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$( centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}
