process FASTME {
    tag "$infile"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastme:2.1.6.1--hec16e2b_1':
        'biocontainers/fastme:2.1.6.1--hec16e2b_1' }"

    input:
    tuple val(meta), path(infile), path(initial_tree)

    output:
    tuple val(meta), path("*.nwk")       , emit: nwk
    tuple val(meta), path("*_stat.txt")  , emit: stats
    tuple val(meta), path("*.matrix.phy"), emit: matrix    , optional: true
    tuple val(meta), path("*.bootstrap") , emit: bootstrap , optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: infile
    def initarg = initial_tree ? "-u $initial_tree" : ''
    def matarg  = task.ext.args =~ "-O" ? "-O ${prefix}.matrix.phy" : ''
    def bootarg = task.ext.args =~ "-B" ? "-B ${prefix}.bootstrap" : ''
    """
    fastme \\
        $args \\
        -i $infile \\
        $initarg \\
        -o ${prefix}.nwk \\
        $matarg \\
        $bootarg \\
        -T $task.cpus


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastme: \$(fastme --version |& sed '1!d ; s/FastME //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: infile
    def mat    = task.ext.args =~ "-O" ? "touch ${prefix}.matrix.phy" : ''
    def boot   = task.ext.args =~ "-B" ? "touch ${prefix}.bootstrap" : ''
    """
    touch ${prefix}.nwk
    touch ${prefix}_stat.txt
    $mat
    $boot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastme: \$(fastme --version |& sed '1!d ; s/FastME //')
    END_VERSIONS
    """
}
