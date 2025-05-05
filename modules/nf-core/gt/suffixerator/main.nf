process GT_SUFFIXERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(fasta)
    val mode

    output:
    tuple val(meta), path("$prefix"), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if ( mode !in [ 'dna', 'protein' ] ) { error "Mode must be one of 'dna', or 'protein'" }
    def args        = task.ext.args     ?: ''
    prefix          = task.ext.prefix   ?: "${meta.id}"
    """
    mkdir \\
        "$prefix"

    gt \\
        suffixerator \\
        "-$mode" \\
        $args \\
        -db $fasta \\
        -indexname "$prefix/suffixerator"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """

    stub:
    if ( mode !in [ 'dna', 'protein' ] ) { error "Mode must be one of 'dna', or 'protein'" }
    def args        = task.ext.args     ?: ''
    prefix          = task.ext.prefix   ?: "${meta.id}"
    def touch_ssp   = mode == "protein" ? "touch $prefix/suffixerator.ssp" : ''
    """
    mkdir \\
        "$prefix"

    touch "$prefix/suffixerator.esq"
    touch "$prefix/suffixerator.prj"
    $touch_ssp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """
}
