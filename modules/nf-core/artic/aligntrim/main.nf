process ARTIC_ALIGNTRIM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/align_trim_samtools:ae17186e78142d18'
        : 'community.wave.seqera.io/library/align_trim_samtools:3883b74edda083a2'}"

    input:
    tuple val(meta), path(samfile), path(scheme_bed), val(normalise_depth)
    val sort_bam

    output:
    tuple val(meta), path("*.primertrimmed*.bam"),    emit: primertrimmed_bam
    tuple val(meta), path("*.align_trim_report.tsv"), emit: align_trim_report
    tuple val(meta), path("*.amp_depth_report.tsv"),  emit: amp_depth_report
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def normalise = normalise_depth ? "--normalise ${normalise_depth}" : ''
    def sort_bam_cmd = sort_bam ? "samtools sort -o ${prefix}.primertrimmed.sorted.bam ${prefix}.primertrimmed.bam && rm ${prefix}.primertrimmed.bam" : ''
    """
    align_trim \\
        --samfile ${samfile} \\
        --output ${prefix}.primertrimmed.bam \\
        --report ${prefix}.align_trim_report.tsv \\
        --amp-depth-report ${prefix}.amp_depth_report.tsv \\
        ${normalise} \\
        ${args} \\
        ${scheme_bed}

    ${sort_bam_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        align_trim: \$(align_trim --version | sed 's/align_trim //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_bam_cmd = sort_bam ? "touch ${prefix}.primertrimmed.sorted.bam" : "touch ${prefix}.primertrimmed.bam"
    """
    echo ${args}
    ${sort_bam_cmd}
    touch ${prefix}.align_trim_report.tsv
    touch ${prefix}.amp_depth_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        align_trim: \$(align_trim --version | sed 's/align_trim //')
    END_VERSIONS
    """
}
