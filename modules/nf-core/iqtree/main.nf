process IQTREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.3.0--h21ec9f0_0' :
        'biocontainers/iqtree:2.3.0--h21ec9f0_0' }"

    input:
    tuple val(meta), path(alignment)
    val constant_sites

    output:
    tuple val(meta), path("*.treefile") , emit: phylogeny
    tuple val(meta), path("*.iqtree")   , emit: report
    tuple val(meta), path("*.mldist")   , emit: mldist, optional: true
    tuple val(meta), path("*.ufboot")   , emit: bootstrap, optional: true
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def fconst_args = constant_sites ? "-fconst $constant_sites" : ''
    def memory      = task.memory.toString().replaceAll(' ', '')
    def prefix      = task.ext.prefix ?: meta.id
    """
    iqtree \\
        $fconst_args \\
        $args \\
        -s $alignment \\
        -pre $prefix \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $memory \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.treefile
    touch ${prefix}.iqtree
    touch ${prefix}.mldist
    touch ${prefix}.ufboot
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """

}
