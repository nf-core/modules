process COBRAMETA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cobra-meta:1.2.3--pyhdfd78af_0':
        'biocontainers/cobra-meta:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(coverage)
    tuple val(meta3), path(query)
    tuple val(meta4), path(bam)
    val assembler
    val mink
    val maxk

    output:
    tuple val(meta), path("${prefix}/COBRA_category_i_self_circular.fasta.gz")                  , emit: self_circular       , optional: true
    tuple val(meta), path("${prefix}/COBRA_category_ii-a_extended_circular_unique.fasta.gz")    , emit: extended_circular   , optional: true
    tuple val(meta), path("${prefix}/COBRA_category_ii-b_extended_partial_unique.fasta.gz")     , emit: extended_partial    , optional: true
    tuple val(meta), path("${prefix}/COBRA_category_ii-c_extended_failed.fasta.gz")             , emit: extended_failed     , optional: true
    tuple val(meta), path("${prefix}/COBRA_category_iii_orphan_end.fasta.gz")                   , emit: orphan_end          , optional: true
    tuple val(meta), path("${prefix}/COBRA_all_assemblies.fasta.gz")                            , emit: all_cobra_assemblies, optional: true
    tuple val(meta), path("${prefix}/COBRA_joining_summary.txt")                                , emit: joining_summary
    tuple val(meta), path("${prefix}/log")                                                      , emit: log
    path "versions.yml"                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cobra-meta \\
        --fasta ${fasta} \\
        --coverage ${coverage} \\
        --query ${query} \\
        --mapping ${bam} \\
        --assembler ${assembler} \\
        --mink ${mink} \\
        --maxk ${maxk} \\
        --threads ${task.cpus} \\
        --output ${prefix} \\
        $args

    gzip ${prefix}/*.fasta
    cat ${prefix}/*fasta.gz > ${prefix}/COBRA_all_assemblies.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    echo "" | gzip > ${prefix}/COBRA_all_assemblies.fasta.gz
    echo "" | gzip > ${prefix}/COBRA_category_i_self_circular.fasta.gz
    echo "" | gzip > ${prefix}/COBRA_category_ii-a_extended_circular_unique.fasta.gz
    echo "" | gzip > ${prefix}/COBRA_category_ii-b_extended_partial_unique.fasta.gz
    echo "" | gzip > ${prefix}/COBRA_category_ii-c_extended_failed.fasta.gz
    echo "" | gzip > ${prefix}/COBRA_category_iii_orphan_end.fasta.gz
    touch ${prefix}/COBRA_joining_summary.txt
    touch ${prefix}/log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """
}
