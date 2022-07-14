process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(database), val(term)

    output:
    tuple val(meta), path("*.xml") , emit: result_xml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    esearch \\
        $args \\
        -db $database \\
        -query $term > | tail -n+3 ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}
