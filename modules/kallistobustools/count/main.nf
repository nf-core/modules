process KALLISTOBUSTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::kb-python=0.26.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.26.3--pyhdfd78af_0' :
        'quay.io/biocontainers/kb-python:0.26.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  t2g
    path  t1c
    path  t2c
    val   workflow_mode
    val   technology

    output:
    tuple val(meta), path ("*.count"), emit: count
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cdna     = t1c ? "-c1 $t1c" : ''
    def introns  = t2c ? "-c2 $t2c" : ''
    """
    kb \\
        count \\
        -t $task.cpus \\
        -i $index \\
        -g $t2g \\
        $cdna \\
        $introns \\
        --workflow $workflow_mode \\
        -x $technology \\
        $args \\
        -o ${prefix}.count \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """
}
