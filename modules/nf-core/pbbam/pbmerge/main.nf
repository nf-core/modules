def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/pbtk/pbmerge

Reason:
This module is no longer fit for purpose because pbbam has been deprecated by PacificBiosciences

"""
process PBBAM_PBMERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbbam:2.4.0--hdcf5f25_1' :
        'biocontainers/pbbam:2.4.0--hdcf5f25_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.pbi"), emit: pbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert true: deprecation_message
    """
    pbmerge \\
        -o ${prefix}.bam \\
        $args \\
        *.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbbam: \$( pbmerge --version | head -n1 | sed 's/pbmerge //' | sed -E 's/ .+//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert true: deprecation_message
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.pbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbbam: \$( pbmerge --version | head -n1 | sed 's/pbmerge //' | sed -E 's/ .+//' )
    END_VERSIONS
    """
}
