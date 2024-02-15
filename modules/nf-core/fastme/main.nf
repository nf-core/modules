process FASTME {
    tag "$infile"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastme:2.1.6.1--hec16e2b_1':
        'biocontainers/fastme:2.1.6.1--hec16e2b_1' }"

    input:
    path infile
    path topo

    output:
    path "*.nwk"        , emit: nwk
    path "*_stat.txt"   , emit: stats
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: infile
    def topoarg = topo ? "-u $topo" : ''
    """
    fastme \\
        $args \\
        -i $infile \\
        $topoarg \\
        -o ${prefix}.nwk \\
        -T $task.cpus


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastme: \$(fastme --version |& sed '1!d ; s/FastME //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: infile
    """
    touch ${prefix}.nwk
    touch ${prefix}_stat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastme: \$(fastme --version |& sed '1!d ; s/FastME //')
    END_VERSIONS
    """
}
