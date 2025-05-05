process ENTREZDIRECT_ESUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(uid), path(uids_file)
    val database

    output:
    tuple val(meta), path("*.xml"), emit: xml
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input = uids_file ? "-input ${uids_file}" : "-id ${uid}"
    if (!uid && !uids_file) error "No input. Valid input: an identifier or a .txt file with identifiers"
    if (uid && uids_file) error "Only one input is required: a single identifier or a .txt file with identifiers"
    """
    esummary \\
        $args \\
        -db $database \\
        $input > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(esummary -version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!uid && !uids_file) error "No input. Valid input: an identifier or a .txt file with identifiers"
    if (uid && uids_file) error "Only one input is required: a single identifier or a .txt file with identifiers"
    """
    touch ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(esummary -version)
    END_VERSIONS
    """
}
