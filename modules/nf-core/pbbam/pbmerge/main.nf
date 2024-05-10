process PBBAM_PBMERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbbam:2.4.0--h8db2425_0' :
        'biocontainers/pbbam:2.4.0--h8db2425_0' }"

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
    """
    touch dummy.bam
    touch dummy.bam.pbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbbam: \$( pbmerge --version | head -n1 | sed 's/pbmerge //' | sed -E 's/ .+//' )
    END_VERSIONS
    """
}
