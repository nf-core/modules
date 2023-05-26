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
    tuple val(meta), path("${prefix}.tsv")                     , emit: coordinates
    tuple val(meta), path("${prefix}.gb*")                     , emit: gbk
    tuple val(meta), path("${prefix}.log")                     , emit: log
    tuple val(meta), path("${prefix}_prophage_information.tsv"), optional:true, emit: information
    tuple val(meta), path("${prefix}_bacteria.fasta")          , optional:true, emit: bacteria_fasta
    tuple val(meta), path("${prefix}_bacteria.gbk")            , optional:true, emit: bacteria_gbk
    tuple val(meta), path("${prefix}_phage.fasta")             , optional:true, emit: phage_fasta
    tuple val(meta), path("${prefix}_phage.gbk")               , optional:true, emit: phage_gbk
    tuple val(meta), path("${prefix}_prophage.gff3")           , optional:true, emit: prophage_gff
    tuple val(meta), path("${prefix}_prophage.tbl")            , optional:true, emit: prophage_tbl
    tuple val(meta), path("${prefix}_prophage.tsv")            , optional:true, emit: prophage_tsv
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Extract GBK file extension, i.e. .gbff, .gbk.gz
    gbk_extension = gbk.getName() - gbk.getSimpleName()

    if ("$gbk" == "${prefix}${gbk_extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    PhiSpy.py \\
        $args \\
        --threads $task.cpus \\
        -p $prefix \\
        -o . \\
        $gbk

    mv ${prefix}_prophage_coordinates.tsv ${prefix}.tsv
    mv ${prefix}_${gbk} ${prefix}${gbk_extension}
    mv ${prefix}_phispy.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    gbk_extension = gbk.getName() - gbk.getSimpleName()

    if ("$gbk" == "${prefix}${gbk_extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.tsv
    touch ${prefix}${gbk_extension}
    touch ${prefix}.log
    touch ${prefix}_prophage_information.tsv
    touch ${prefix}_bacteria.fasta
    touch ${prefix}_bacteria.gbk
    touch ${prefix}_phage.fasta
    touch ${prefix}_phage.gbk
    touch ${prefix}_prophage.gff3
    touch ${prefix}_prophage.tbl
    touch ${prefix}_prophage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """
}
