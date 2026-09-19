process NANOMONSV_PARSE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanomonsv:0.8.0--pyhdfd78af_0':
        'quay.io/biocontainers/nanomonsv:0.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.insertion.sorted.bed.gz")          , emit: insertions
    tuple val(meta), path("${prefix}.insertion.sorted.bed.gz.tbi")      , emit: insertions_index
    tuple val(meta), path("${prefix}.deletion.sorted.bed.gz")           , emit: deletions
    tuple val(meta), path("${prefix}.deletion.sorted.bed.gz.tbi")       , emit: deletions_index
    tuple val(meta), path("${prefix}.rearrangement.sorted.bedpe.gz")    , emit: rearrangements
    tuple val(meta), path("${prefix}.rearrangement.sorted.bedpe.gz.tbi"), emit: rearrangements_index
    tuple val(meta), path("${prefix}.bp_info.sorted.bed.gz")            , emit: bp_info
    tuple val(meta), path("${prefix}.bp_info.sorted.bed.gz.tbi")        , emit: bp_info_index
    tuple val("${task.process}"), val('nanomonsv'), eval("nanomonsv --version 2>&1 | sed 's/^nanomonsv //'") , topic: versions, emit: versions_nanomonsv
    tuple val("${task.process}"), val('mafft'), eval("mafft --version 2>&1 | sed 's/^v//; s/ (.*//'")        , topic: versions, emit: versions_mafft
    tuple val("${task.process}"), val('racon'), eval("racon --version 2>&1 | sed 's/^v//'")                  , topic: versions, emit: versions_racon
    tuple val("${task.process}"), val('tabix'), eval("tabix --version 2>&1 | sed '1!d;s/^tabix (htslib) //'"), topic: versions, emit: versions_tabix
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version 2>&1 | sed '1!d;s/^bgzip (htslib) //'"), topic: versions, emit: versions_bgzip
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //g'")              , topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args   ?: ''

    """
    nanomonsv parse \\
        ${args} \\
        ${bam} \\
        ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.insertion.sorted.bed.gz
    touch ${prefix}.insertion.sorted.bed.gz.tbi
    echo "" | gzip > ${prefix}.deletion.sorted.bed.gz
    touch ${prefix}.deletion.sorted.bed.gz.tbi
    echo "" | gzip > ${prefix}.rearrangement.sorted.bedpe.gz
    touch ${prefix}.rearrangement.sorted.bedpe.gz.tbi
    echo "" | gzip > ${prefix}.bp_info.sorted.bed.gz
    touch ${prefix}.bp_info.sorted.bed.gz.tbi
    """
}
