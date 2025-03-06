process CIRCULARMAPPER_REALIGNSAMFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circularmapper:1.93.5--h4a94de4_1':
        'biocontainers/circularmapper:1.93.5--h4a94de4_1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), val(elongation_factor)
    tuple val(meta4), path(elongated_chr_list)
    // NOTE: The elongated_chr_list is not used in the script, but is an implicit input that realignsamfile requires when using the `-f true` option.
    //          In its absence, when `-f true` is set, realignsamfile will remove all @SQ tags from the BAM header, breaking the bamfile.

    output:
    tuple val(meta), path("*_realigned.bam")    , emit: bam
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.93.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    realignsamfile \\
        -Xmx${task.memory.toGiga()}g \\
        ${args} \\
        -e ${elongation_factor} \\
        -i ${bam} \\
        -r ${fasta}

    ## realignsamfile has a hardcoded output name. Rename if necessary to use prefix.
    if [[ "${bam.getBaseName()}_realigned.bam" != "${prefix}_realigned.bam" ]]; then
        mv ${bam.getBaseName()}_realigned.bam ${prefix}_realigned.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CircularMapper: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.93.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_realigned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CircularMapper: ${VERSION}
    END_VERSIONS
    """
}
