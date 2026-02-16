process HOMER_FINDPEAKS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // singularity build url: https://wave.seqera.io/view/builds/bd-9c603739ae7d4fd3_1
    // docker build url: https://wave.seqera.io/view/builds/bd-08c7bb832e96c6bd_1
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:9c603739ae7d4fd3'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:08c7bb832e96c6bd'}"

    input:
    val(style)
    tuple val(meta), path(tagDir)
    tuple val(meta2), path(controlTagDir)
    path uniqmap

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.gtf"), optional: true, emit: gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def uniqmap_flag = uniqmap ? "-uniqmap ${uniqmap}" : ""
    def control_flag = controlTagDir ? "-i ${controlTagDir}" : ""
    def control_suffix = controlTagDir ? "_vs_${meta2.id}" : ""
    def VERSION = '5.1'

    // Validate style
    def valid_styles = ['factor', 'histone', 'groseq', 'tss', 'dnase', 'super', 'mC']
    if (!valid_styles.contains(style)) {
        error "style must be one of: ${valid_styles.join(', ')}"
    }

    // Determine output suffix based on style
    def output_suffix_map = [
        'factor': 'peaks.txt',
        'histone': 'regions.txt',
        'groseq': 'transcripts.txt',
        'tss': 'tss.txt',
        'dnase': 'peaks.txt',
        'super': 'superEnhancers.txt',
        'mC': 'regions.txt'
    ]
    def output_suffix = output_suffix_map[style]
    def output_filename = "${prefix}${control_suffix}_${output_suffix}"

    def groseq_gtf_arg = style == "groseq" ? "-gtf ${prefix}${control_suffix}_transcripts.gtf" : ""

    """
    findPeaks \\
        ${tagDir} \\
        -style ${style} \\
        ${control_flag} \\
        ${args} \\
        ${groseq_gtf_arg} \\
        -o ${output_filename} \\
        ${uniqmap_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control_suffix = controlTagDir ? "_vs_${meta2.id}" : ""
    def VERSION = '5.1'

    def output_suffix_map = [
        'factor': 'peaks.txt',
        'histone': 'regions.txt',
        'groseq': 'transcripts.txt',
        'tss': 'tss.txt',
        'dnase': 'peaks.txt',
        'super': 'superEnhancers.txt',
        'mC': 'regions.txt'
    ]
    def output_suffix = output_suffix_map[style]
    def output_filename = "${prefix}${control_suffix}_${output_suffix}"
    def gtf_filename = "${prefix}${control_suffix}_transcripts.gtf"

    """
    touch ${output_filename}
    ${style == 'groseq' ? "touch ${gtf_filename}" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """
}
