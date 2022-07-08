process ENTREZDIRECT_ESUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta)
    val uid
    path ids
    val database

    output:
    tuple val(meta), path("*.xml"), emit: xml_esummary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def std_input = uid ? "-id ${uid}" : null
    def file_input = ids ? "-input ${ids}" : null
    // std_input: single uid
    if( std_input != null && file_input == null)
    """
    esummary \\
        $args \\
        -db $database \\
        $std_input | cat | grep '<' > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(echo \$(esummary --help | head -1 | cut -d' ' -f2 2>&1) | sed 's/^esummary //; s/Using.*\$//' ))
    END_VERSIONS
    """
    // file_input: list of ids, one id per line
    else if( file_input != null &&  std_input == null)
    """
    esummary \\
        $args \\
        -db $database \\
        $file_input | cat | grep '<' > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esummary: \$(echo \$(esummary --help | head -1 | cut -d' ' -f2 2>&1) | sed 's/^esummary //; s/Using.*\$//' ))
    END_VERSIONS
    """
    else
       error "Invalid input: provide a single unique identifier or a file of identifiers."
}
