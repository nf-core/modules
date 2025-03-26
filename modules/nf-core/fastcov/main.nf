process FASTCOV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ? 'community.wave.seqera.io/library/fastcov:0.1.3--84def91a6ef27f61': 'oras://community.wave.seqera.io/library/fastcov:0.1.3--e8503352127532b6' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("${meta.id}.${task.ext.file_ext ?: 'png'}"), emit: coverage_plot
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_ext = task.ext.file_ext ?: 'png'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    fastcov.py -o ${prefix}.${file_ext} ${args} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcov: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def file_ext = task.ext.file_ext ?: 'png'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.${file_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcov: $VERSION
    END_VERSIONS
    """
}
