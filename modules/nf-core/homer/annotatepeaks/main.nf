process HOMER_ANNOTATEPEAKS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // singularity build url: https://wave.seqera.io/view/builds/bd-9c603739ae7d4fd3_1
    // docker build url: https://wave.seqera.io/view/builds/bd-08c7bb832e96c6bd_1
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:9c603739ae7d4fd3'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:08c7bb832e96c6bd'}"

    input:
    tuple val(meta), path(peak)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*annotatePeaks.txt"), emit: txt
    tuple val(meta), path("*annStats.txt"), emit: stats, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    annotatePeaks.pl \\
        ${peak} \\
        ${fasta} \\
        ${args} \\
        -gtf ${gtf} \\
        -cpu ${task.cpus} \\
        > ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """
}
