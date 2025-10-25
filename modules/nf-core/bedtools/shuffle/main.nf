process BEDTOOLS_SHUFFLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3':
        'biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(meta) , path(intervals)
    tuple val(meta2), path(chrom_sizes)
    path exclude_file
    path include_file

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def exclude_arg = exclude_file ? "-excl $exclude_file" : ''
    def include_arg = include_file ? "-incl $include_file" : ''
    if ("$intervals" == "${prefix}.bed") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    bedtools \\
        shuffle \\
        $args \\
        -i $intervals \\
        -g $chrom_sizes \\
        $exclude_arg \\
        $include_arg \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed '1!d ; s/bedtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$intervals" == "${prefix}.bed") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed '1!d ; s/bedtools //')
    END_VERSIONS
    """
}
