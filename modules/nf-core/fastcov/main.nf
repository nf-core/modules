process FASTCOV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ? 
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe768c866b4ae4f2c8948ea7a274ce66d590eedb1cf967495dbd0fb84643a7e2/data': 
        'community.wave.seqera.io/library/fastcov:0.1.3--84def91a6ef27f61' }"

    input:
    tuple val(meta), path(bam), path(index)
    val(file_ext)

    output:
    tuple val(meta), path("${prefix}.${file_ext}"), emit: coverage_plot
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = file_ext.isEmpty() ? 'png' : file_ext
    def VERSION = '0.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    fastcov.py \\
      -o ${prefix}.${file_ext} \\
      ${args} \\
      ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcov: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = file_ext.isEmpty() ? 'png' : file_ext
    def VERSION = '0.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.${file_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcov: $VERSION
    END_VERSIONS
    """
}
