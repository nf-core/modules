process ENTREZDIRECT_ESUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::entrez-direct=13.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:13.9--pl526h375a9b1_1':
        'quay.io/biocontainers/entrez-direct:13.9--pl526h375a9b1_1' }"

    input:
    tuple val(meta), val(uid)
    val database

    output:
    tuple val(meta), path("*.xml"), emit: xml_esummary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    esummary \\
        -db $database \\
        -id $uid \\
        $args \\
        > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(echo \$(esummary --help | head -1 | cut -d' ' -f2 2>&1) | sed 's/^esummary //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
