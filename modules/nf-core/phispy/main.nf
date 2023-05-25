process PHISPY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::phispy=4.2.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phispy:4.2.21--py310h30d9df9_1':
        'biocontainers/phispy:4.2.21--py310h30d9df9_1' }"

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: coordinates
    tuple val(meta), path("${prefix}.gb*"), emit: gbk
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Extract GBK file extension, i.e. .gbff, .gbk.gz
    def gbk_extension = gbk.getName() - gbk.getSimpleName()

    if ("$gbk" == "${prefix}${gbk_extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    PhiSpy.py \\
        $args \\
        --threads $task.cpus \\
        -p $prefix \\
        -o $prefix \\
        $gbk

    mv ${prefix}/${prefix}_prophage_coordinates.tsv ${prefix}.tsv
    mv ${prefix}/${prefix}_${gbk} ${prefix}${gbk_extension}
    mv ${prefix}/${prefix}_phispy.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.gbk
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """
}
