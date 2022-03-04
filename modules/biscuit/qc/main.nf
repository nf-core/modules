process BISCUIT_QC {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.0.2.20220113--h81a5ba2_0':
        'quay.io/biocontainers/biscuit:1.0.2.20220113--h81a5ba2_0' }"

    input:
    tuple val(meta), path(bam)
    path(index)

    output:
    tuple val(meta), path("*.txt"), emit: biscuit_qc_reports
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def se = meta.single_end ? "-s" : ""

    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

    biscuit qc \\
        $se \\
        \$INDEX \\
        $bam \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$(echo \$(biscuit version 2>&1) | sed 's/^.*BISCUIT Version: //; s/Using.*\$//')
    END_VERSIONS
    """
}
