process FASTANI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.32--he1c1bb9_0' :
        'biocontainers/fastani:1.32--he1c1bb9_0' }"

    input:
    tuple val(meta), path(query)
    tuple val(meta2), path(reference)
    path ql
    path rl

    output:
    tuple val(meta), path("*.txt")  , emit: ani
    tuple val(meta), path("*.visual")   , optional:true, emit: visual
    tuple val(meta), path("*.matrix")   , optional:true, emit: matrix
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix2 = task.ext.prefix2 ?: "${meta2.id}"
    def input_query = query ? "-q ${query}": "--ql ${ql}"
    def input_reference = reference ? "-r ${reference}": "--rl ${rl}"
    def out = query ?: ( reference ? "-o ${prefix}_v_${prefix2}": "-o ${prefix}_v_all"): ( reference ? "-o all_v_${prefix2}": "-o all_v_all" )

    """
    fastANI \\
        $input_query \\
        $input_reference \\
        --threads $task.cpus \\
        $out
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix2 = task.ext.prefix2 ?: "${meta2.id}"
    def input_query = query ? "-q ${query}": "--ql ${ql}"
    def input_reference = reference ? "-r ${reference}": "--rl ${rl}"
    def out = query ?: ( reference ? "-o ${prefix}_v_${prefix2}": "-o ${prefix}_v_all"): ( reference ? "-o all_v_${prefix2}": "-o all_v_all" )

    """
    touch ${out}.visual
    touch ${out}.txt
    touch ${out}.matrix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}