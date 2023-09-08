process NANOMONSV_PARSE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::nanomonsv=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanomonsv:0.5.0--pyhdfd78af_0':
        'biocontainers/nanomonsv:0.5.0--pyhdfd78af_0' }"

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
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    nanomonsv parse ${args} ${bam} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(echo \$(nanomonsv --version 2>&1) | sed 's/^nanomonsv //')
        mafft: \$(echo \$(mafft --version 2>&1) | sed 's/^v//; s/ (.*//')
        racon: \$(echo \$(racon --version 2>&1) | sed 's/^v//')
        tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^tabix (htslib) //; s/ Copyright.*//')
        bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^bgzip (htslib) //; s/ Copyright.*//')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
