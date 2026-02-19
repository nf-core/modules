process FASTANI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.34--hb66fcc3_7' :
        'biocontainers/fastani:1.34--hb66fcc3_7' }"

    input:
    tuple val(meta),  path(query)
    tuple val(meta2), path(reference)
    path(ql)
    path(rl)

    output:
    tuple val(meta), path("*.txt")   , emit: ani
    tuple val(meta), path("*.visual"), optional:true, emit: visual
    tuple val(meta), path("*.matrix"), optional:true, emit: matrix
    tuple val("${task.process}"), val("fastani"), eval('fastANI --version 2>&1 | head -1 | sed "s/version\\ //"'), topic: versions, emit: versions_fastani

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix ?: ( meta.id ?: 'all' )
    def input_query     = query ? "-q ${query}": "--ql ${ql}"
    def input_reference = reference ? "-r ${reference}": "--rl ${rl}"
    """
    fastANI \\
        $input_query \\
        $input_reference \\
        --threads $task.cpus \\
        -o ${prefix}.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix  ?: ( meta.id  ?: 'all')
    """
    echo $args

    touch ${prefix}.visual
    touch ${prefix}.txt
    touch ${prefix}.matrix
    """
}
