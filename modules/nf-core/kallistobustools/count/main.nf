process KALLISTOBUSTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.29.1--pyhdfd78af_0' :
        'biocontainers/kb-python:0.29.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  t2g
    path  t1c
    path  t2c
    path  whitelist
    val   technology
    val   workflow_mode

    output:
    tuple val(meta), path ("*.count")   , emit: count
    path "versions.yml"                 , emit: versions
    path "*.count/*/*.mtx"              , emit: matrix //Ensure that kallisto finished and produced outputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def cdna      = t1c ? "-c1 $t1c" : ''
    def intron    = t2c ? "-c2 $t2c" : ''
    def whitelist = whitelist ? "-w $whitelist" : ''
    def memory    = task.memory.toGiga() - 1
    """
    kb \\
        count \\
        -t $task.cpus \\
        -i $index \\
        -g $t2g \\
        $cdna \\
        $intron \\
        $whitelist \\
        -x $technology \\
        --workflow $workflow_mode \\
        $args \\
        -o ${prefix}.count \\
        -m ${memory}G \\
        ${reads.join( " " )}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.count/counts_unfiltered/
    touch ${prefix}.count/counts_unfiltered/cells_x_genes.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """
}
