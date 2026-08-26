process PHISPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phispy:4.2.21--py310h30d9df9_1':
        'quay.io/biocontainers/phispy:4.2.21--py310h30d9df9_1' }"

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${prefix}.tsv")                     , emit: coordinates
    tuple val(meta), path("${prefix}.gb*")                     , emit: gbk
    tuple val(meta), path("${prefix}.log")                     , emit: log
    tuple val(meta), path("${prefix}_prophage_information.tsv"), emit: information   , optional:true
    tuple val(meta), path("${prefix}_bacteria.fasta")          , emit: bacteria_fasta, optional:true
    tuple val(meta), path("${prefix}_bacteria.gbk")            , emit: bacteria_gbk  , optional:true
    tuple val(meta), path("${prefix}_phage.fasta")             , emit: phage_fasta   , optional:true
    tuple val(meta), path("${prefix}_phage.gbk")               , emit: phage_gbk     , optional:true
    tuple val(meta), path("${prefix}_prophage.gff3")           , emit: prophage_gff  , optional:true
    tuple val(meta), path("${prefix}_prophage.tbl")            , emit: prophage_tbl  , optional:true
    tuple val(meta), path("${prefix}_prophage.tsv")            , emit: prophage_tsv  , optional:true
    tuple val("${task.process}"), val('phispy'), eval('PhiSpy.py --version 2>&1'), topic: versions, emit: versions_phispy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Extract GBK file extension, i.e. .gbff, .gbk.gz
    gbk_extension = gbk.getName() - gbk.getSimpleName()

    if ("${gbk}" == "${prefix}${gbk_extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    PhiSpy.py \\
        ${args} \\
        --threads ${task.cpus} \\
        -p ${prefix} \\
        -o . \\
        ${gbk}

    mv ${prefix}_prophage_coordinates.tsv ${prefix}.tsv
    mv ${prefix}_${gbk} ${prefix}${gbk_extension}
    mv ${prefix}_phispy.log ${prefix}.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    gbk_extension = gbk.getName() - gbk.getSimpleName()
    gbl_create_cmd = gbk_extension.endsWith(".gz") ? 'echo "" | gzip >' : "touch"

    if ("${gbk}" == "${prefix}${gbk_extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.tsv
    ${gbl_create_cmd} ${prefix}${gbk_extension}
    touch ${prefix}.log
    touch ${prefix}_prophage_information.tsv
    touch ${prefix}_bacteria.fasta
    touch ${prefix}_bacteria.gbk
    touch ${prefix}_phage.fasta
    touch ${prefix}_phage.gbk
    touch ${prefix}_prophage.gff3
    touch ${prefix}_prophage.tbl
    touch ${prefix}_prophage.tsv
    """
}
