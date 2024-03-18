process IQTREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.1.4_beta--hdcc8f71_0' :
        'biocontainers/iqtree:2.1.4_beta--hdcc8f71_0' }"

    input:
    tuple val(meta), path(alignment)
    val constant_sites

    output:
    path "*.treefile"   , emit: phylogeny
    path "versions.yml" , emit: versions

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
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $memory \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
