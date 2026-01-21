process LOCALCDSEARCH_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/local-cd-search:0.3.0--pyhdfd78af_0' :
        'biocontainers/local-cd-search:0.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    val sites

    output:
    tuple val(meta), path("*_results.tsv"), emit: result
    tuple val(meta), path("*_sites.tsv")  , emit: annot_sites, optional: true
    tuple val("${task.process}"), val('local-cd-search'), eval("echo ${VERSION}"), topic: versions, emit: versions_localcdsearch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    VERSION = '0.3.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def prefix = task.ext.prefix ?: "${meta.id}"
    def val_flag = sites ? "--sites-output ${prefix}_sites.tsv" : ''
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''
    """
    $uncompress_input

    local-cd-search \\
        annotate \\
        $args \\
        $val_flag \\
        --threads $task.cpus \\
        ${fasta_name} \\
        ${prefix}_results.tsv \\
        ${db}
    """

    stub:
    def args = task.ext.args ?: ''
    VERSION = '0.3.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}_results.tsv
    ${sites ? "touch ${prefix}_sites.tsv" : ''}
    """
}
