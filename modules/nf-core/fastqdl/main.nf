process FASTQDL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6a69a81750a1808471528b20358f5f03b555d0333d622e2f3e47e33700a6718/data':
        'community.wave.seqera.io/library/fastq-dl:4.0.1--fd5420f75c055a44' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fastq.gz")       , emit: fastq
    tuple val(meta), path("*-run-info.tsv")   , emit: runinfo
    tuple val(meta), path("*-run-mergers.tsv"), emit: runmergers, optional: true
    tuple val("${task.process}"), val('fastq-dl'), eval('fastq-dl --version |& sed "s/.* //"'), emit: versions_fastqdl, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq-dl \\
        ${args} \\
        --prefix ${prefix} \\
        --accession ${accession} \\
        --cpus ${task.cpus} \\
        --outdir .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${accession}.fastq.gz
    echo "" | gzip > ${accession}_1.fastq.gz
    echo "" | gzip > ${accession}_2.fastq.gz
    touch ${prefix}-run-info.tsv
    """
}
