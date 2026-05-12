process GT_SUFFIXERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'quay.io/biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(fasta)
    val mode

    output:
    tuple val(meta), path("$prefix"), emit: index
    tuple val("${task.process}"), val('genometools'), eval("gt --version | sed '1!d;s/gt (GenomeTools) //'"), emit: versions_gt, topic: versions

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
        -${mode} \\
        ${args} \\
        -db ${fasta} \\
        -indexname ${prefix}/suffixerator
    """

    stub:
    if ( mode !in [ 'dna', 'protein' ] ) { error "Mode must be one of 'dna', or 'protein'" }
    prefix          = task.ext.prefix   ?: "${meta.id}"
    def touch_ssp   = mode == "protein" ? "touch ${prefix}/suffixerator.ssp" : ''
    """
    mkdir ${prefix}

    touch ${prefix}/suffixerator.esq
    touch ${prefix}/suffixerator.prj
    ${touch_ssp}
    """
}
