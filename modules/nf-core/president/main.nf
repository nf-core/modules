process PRESIDENT {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/president:0.6.8--pyhdfd78af_0' :
        'biocontainers/president:0.6.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(reference)
    val compress

    output:
    tuple val(meta), path("${prefix}_valid.fasta*")  , emit: valid_fasta
    tuple val(meta), path("${prefix}_invalid.fasta*"), emit: invalid_fasta
    tuple val(meta), path("*.tsv")                   , emit: report
    tuple val(meta), path("*.log")                   , emit: log
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${fasta}" == "${reference}") error "Input and reference names are the same!"
    if ("${fasta}" == "${prefix}_valid.fasta") error "Input and output file names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("${fasta}" == "${prefix}_invalid.fasta") error "Input and output file names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("${fasta}" == "${prefix}_valid.fasta.gz" && compress) error "Input and output file names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("${fasta}" == "${prefix}_invalid.fasta.gz" && compress) error "Input and output file names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    president \\
        --query $fasta \\
        --reference $reference \\
        --path . \\
        --threads $task.cpus \\
        --prefix "${prefix}_" \\
        $args

    if [ "$compress" = true ] ; then
        gzip ${prefix}*.fasta;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        president: \$(president --version |& sed '1!d ; s/president v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    compress_command = compress ? "gzip ${prefix}*.fasta" : ""

    """
    touch ${prefix}_report.tsv
    touch ${prefix}_president_logger.log
    touch ${prefix}_valid.fasta
    touch ${prefix}_invalid.fasta
    $compress_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        president: \$(president --version |& sed '1!d ; s/president v//')
    END_VERSIONS
    """
}
