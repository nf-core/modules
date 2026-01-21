process ENTREZDIRECT_XTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
    'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), path(xml_input)
    val pattern
    val element
    val sep

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $xml_input | xtract \\
        -pattern $pattern \\
        -tab $sep \\
        -element $element \\
        $args \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xtract: \$(xtract -version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xtract: \$(xtract -version)
    END_VERSIONS
    """
}
