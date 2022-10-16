process METAPHLAN3_MERGEMETAPHLANTABLES {
    label 'process_single'

    conda (params.enable_conda ? 'bioconda::metaphlan=3.0.12' : null)
    def container_image = "metaphlan:3.0.12--pyhb7b1952_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    input:
    tuple val(meta), path(profiles)

    output:
    tuple val(meta), path("${prefix}.txt") , emit: txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_metaphlan_tables.py \\
        $args \\
        -o ${prefix}.txt \\
        ${profiles}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan3: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
