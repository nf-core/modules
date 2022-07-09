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
    def input_id = uid ? "-id ${uid}" : ''
    def input_file = uids_file.name != 'NO_FILE' ? "-input $uids_file" : ''
    if(input_id && !input_file)
    """
    esummary \\
        $args \\
        -db $database \\
        $input_id | cat | grep '<' > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(echo \$(esummary --help | head -1 | cut -d' ' -f2 2>&1) | sed 's/^esummary //; s/Using.*\$//' ))
    END_VERSIONS
    """
    else if(input_file && !input_id)
    """
    esummary \\
        $args \\
        -db $database \\
        $input_file | cat | grep '<' > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(echo \$(esummary --help | head -1 | cut -d' ' -f2 2>&1) | sed 's/^esummary //; s/Using.*\$//' ))
    END_VERSIONS
    """
    else
        error "Invalid input: provide an unique identifier or a file with identifiers."
}
