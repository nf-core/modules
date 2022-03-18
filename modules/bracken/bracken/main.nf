process BRACKEN_BRACKEN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bracken=2.6.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bracken:2.6.2--py39hc16433a_0':
        'quay.io/biocontainers/bracken:2.6.2--py39hc16433a_0' }"

    input:
    tuple val(meta), path(kraken_report), path(database), val(taxonomic_level)

    output:
    tuple val(meta), path(report), emit: reports
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-t 10'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bracken_version = '2.6.1'
    report = "${prefix}_${taxonomic_level}.tsv"
    """
    bracken \\
        ${args} \\
        -d "${database}" \\
        -i "${kraken_report}" \\
        -o "${report}" \\
        -r ${meta.read_length} \\
        -l ${taxonomic_level}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${bracken_version}
    END_VERSIONS
    """
}
