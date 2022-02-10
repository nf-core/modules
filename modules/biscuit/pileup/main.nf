process BISCUIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113 bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-44134f7dad96451eeaecf3505064bb0a4ff131aa:16678cc03bf715c528c7d2db1f0fb570d2358ce8-0':
        'quay.io/biocontainers/mulled-v2-44134f7dad96451eeaecf3505064bb0a4ff131aa:16678cc03bf715c528c7d2db1f0fb570d2358ce8-0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    biscuit pileup \\
        $args \\

        $index \\
        | bgzip $args2 > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
