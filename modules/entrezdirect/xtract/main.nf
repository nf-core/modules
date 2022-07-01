process ENTREZDIRECT_XTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::entrez-direct=13.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/entrez-direct:13.9--pl526h375a9b1_1':
    'quay.io/biocontainers/entrez-direct:13.9--pl526h375a9b1_1' }"

    input:
    tuple val(meta), path(xml_input)
    val pattern
    val element

    output:
    tuple val(meta), path("*.xtract.table.txt"), emit: xtract_table
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $xml_input | xtract \\
        -pattern $pattern \\
        -element $element \\
        $args \\
        > ${prefix}.xtract.table.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        entrezdirect: \$(echo \$(xtract -help | head -1 | cut -d' ' -f2 2>&1) | sed 's/^xtract //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
