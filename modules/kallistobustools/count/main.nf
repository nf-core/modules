process KALLISTOBUSTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::kb-python=0.27.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.27.2--pyhdfd78af_0' :
        'quay.io/biocontainers/kb-python:0.27.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  t2g
    path  t1c
    path  t2c
    val   technology

    output:
    tuple val(meta), path ("*.count"), emit: count
    path "versions.yml"              , emit: versions
    path "*.count/*/*.mtx"           , emit: matrix //Ensure that kallisto finished and produced outputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def cdna    = t1c ? "-c1 $t1c" : ''
    def introns = t2c ? "-c2 $t2c" : ''
    def memory  = task.memory.toGiga() - 1
    """
    kb \\
        count \\
        -t $task.cpus \\
        -i $index \\
        -g $t2g \\
        $cdna \\
        $introns \\
        -x $technology \\
        $args \\
        -o ${prefix}.count \\
        ${reads.join( " " )} \\
        -m ${memory}G

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """
}
