process WHAMG {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::wham=1.8.0 bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fb3127231a8f88daa68d5eb0eec49593cd98b440:e6a1c182ebdfe372b5e3cdc05a6cece64cef7274-0':
        'biocontainers/mulled-v2-fb3127231a8f88daa68d5eb0eec49593cd98b440:e6a1c182ebdfe372b5e3cdc05a6cece64cef7274-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , emit: tbi
    tuple val(meta), path("*.txt")          , emit: graph, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    whamg \\
        -x ${task.cpus} \\
        -f ${bam} \\
        -a ${fasta} \\
        ${args} \\
        | bgzip ${args2} --threads ${task.cpus} --stdout > ${prefix}.vcf.gz

        tabix ${args3} ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whamg: \$(echo \$(whamg 2>&1 | grep Version | sed 's/^Version: v//; s/-.*\$//' ))
    END_VERSIONS
    """
}
