process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    def container_image = "/entrez-direct:16.2--he881be0_1"
                                                       container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), val(term)
    val database

    output:
    tuple val(meta), path("*.xml") , emit: xml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    esearch \\
        -db $database \\
        -query $term \\
        $args > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}
