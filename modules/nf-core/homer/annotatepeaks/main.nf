process HOMER_ANNOTATEPEAKS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0fe4a3875b78dce3c66b43fb96489769cc32e55e329e2525d2af09096af2252a/data'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:a8f4c58755bb281b'}"

    input:
    tuple val(meta), path(peak)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*annotatePeaks.txt"), emit: txt
    tuple val(meta), path("*annStats.txt"), emit: stats, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('homer'), val("4.11"), emit: versions_homer, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotatePeaks.pl \\
        ${peak} \\
        ${fasta} \\
        ${args} \\
        -gtf ${gtf} \\
        -cpu ${task.cpus} \\
        > ${prefix}.annotatePeaks.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotatePeaks.txt
    """
}
