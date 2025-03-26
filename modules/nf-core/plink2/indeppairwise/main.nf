process PLINK2_INDEPPAIRWISE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0':
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(plink_file1), path(plink_file2), path(plink_file3)
    val(win)
    val(step)
    val(r2)

    output:
    tuple val(meta), path("*.prune.in")  , emit: prune_in
    tuple val(meta), path("*.prune.out") , emit: prune_out
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = ([plink_file1.extension, plink_file2.extension, plink_file3.extension].any { it.contains('pgen') }) ? '--pfile' : '--bfile'
    def input = "${plink_file1.getBaseName()}"
    """
    plink2 \\
        $mode $input \\
        $args \\
        --indep-pairwise $win $step $r2 \\
        --threads $task.cpus \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    touch ${prefix}.prune.in
    touch ${prefix}.prune.out
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version)
    END_VERSIONS
    """
}
