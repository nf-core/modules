process PBBAM_PBMERGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pbbam=1.7.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbbam:1.7.0--h058f120_1' :
        'quay.io/biocontainers/pbbam:1.7.0--h058f120_1' }"

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
        pbbam: \$( pbmerge --version|sed 's/pbmerge //' )
    END_VERSIONS
    """
}
