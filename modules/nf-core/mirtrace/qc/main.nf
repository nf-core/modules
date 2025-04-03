process MIRTRACE_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirtrace:1.0.1--0':
        'biocontainers/mirtrace:1.0.1--0' }"

    input:
    tuple val(meta), path(reads), path(mirtrace_config)
    val(mirtrace_species)

    output:
    tuple val(meta), path ("*.html")                                                 , emit: html
    tuple val(meta), path ("*.json")                                                 , emit: json
    tuple val(meta), path ("*.tsv")                                                  , emit: tsv
    tuple val(meta), path ("qc_passed_reads.all.collapsed/*.{fa,fasta}")             , emit: all_fa
    tuple val(meta), path ("qc_passed_reads.rnatype_unknown.collapsed/*.{fa,fasta}") , emit: rnatype_unknown_fa
    path "versions.yml"                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mirtrace_mode = mirtrace_config ? "--config ${mirtrace_config}": "${reads}"

    """
    mirtrace qc  \\
        --species ${mirtrace_species} \\
        --write-fasta \\
        --output-dir . \\
        --force \\
        ${mirtrace_mode}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtrace: \$(echo \$(mirtrace -v))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.html
    touch ${prefix}.json
    touch ${prefix}.tsv

    mkdir -p qc_passed_reads.all.collapsed
    mkdir -p qc_passed_reads.rnatype_unknown.collapsed

    touch qc_passed_reads.all.collapsed/${prefix}.fa
    touch qc_passed_reads.rnatype_unknown.collapsed/${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtrace: \$(echo \$(mirtrace -v))
    END_VERSIONS
    """
}
