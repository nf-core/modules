process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::entrez-direct=16.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

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
