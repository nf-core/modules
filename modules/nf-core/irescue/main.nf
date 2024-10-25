process IRESCUE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/irescue:1.1.2--pyhdfd78af_0':
        'biocontainers/irescue:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    val genome
    path bed

    output:
    tuple val(meta), path("${prefix}")            , emit: results
    tuple val(meta), path("${prefix}/counts")     , emit: counts
    tuple val(meta), path("${prefix}/irescue.log"), emit: log
    tuple val(meta), path("${prefix}/tmp")        , emit: tmp, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reference = bed ? "--regions $bed" : ''
    def genome_assembly = reference ? '' : "--genome $genome"
    """
    mkdir -p $prefix

    irescue \\
        --bam $bam \\
        $reference \\
        $genome_assembly \\
        --outdir $prefix \\
        --threads $task.cpus \\
        $args 2> >(tee -a ${prefix}/irescue.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irescue: \$(irescue --version |& sed '1!d ; s/IRescue //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/counts
    touch \\
        ${prefix}/counts/matrix.mtx \\
        ${prefix}/counts/barcodes.tsv \\
        ${prefix}/counts/features.tsv \\
        ${prefix}/irescue.log
    gzip ${prefix}/counts/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irescue: \$(irescue --version |& sed '1!d ; s/IRescue //')
    END_VERSIONS
    """
}
