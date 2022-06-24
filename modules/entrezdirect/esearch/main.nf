process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::entrez-direct=13.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:13.9--pl526h375a9b1_1':
        'quay.io/biocontainers/entrez-direct:13.9--pl526h375a9b1_1' }"

    input:
    tuple val(meta), val(database), val(term)
    val sort_by
    val date_type
    val min_date
    val max_date
    val spell_check

    output:
    tuple val(meta), path('*.esearch.xml.txt') , emit: xml_output
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_by = meta.sort_by ? "-sort ${meta.sort_by}" : ''
    def date_type = meta.date_type ? "-datetype ${meta.date_type}" : ''
    def min_date = meta.min_date ? "-mindate ${meta.min_date}" : ''
    def max_date = meta.max_date ? "-maxdate ${meta.max_date}" : ''
    def spell_check = meta.spell_check ? "-spell" : ''
    def args = task.ext.args ?: "${sort_by} ${date_type} ${min_date} ${max_date} ${spell_check}"
    """
    esearch \\
        -db $database \\
        -query $term \\
        ${args} \\
        > ${prefix}.esearch.xml.txt 2> ${prefix}.esearch.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        entrezdirect: \$(echo \$(esearch -version 2>&1) | sed 's/^esearch //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
