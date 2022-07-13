process ENTREZDIRECT_ESUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(uid), path(uids_file)
    val database

    output:
    tuple val(meta), path("*.xml"), emit: xml_esummary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input = uids_file ? "-input ${uids_file}" : "-id ${uid}"
    if (!uid && !uids_file) error "No input. Valid input: an identifier or a .txt file with identifiers"
    if (uid && uids_file) error "Only one input is required: a single identifier or a .txt file with identifiers"
    // use of grep is to ensure a clean XML file. Otherwise an irrelevant Perl compilation error ends up in the XML
    """
    esummary \\
        $args \\
        -db $database \\
        $input | grep '<' > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(esummary -version)
    END_VERSIONS
    """
}
