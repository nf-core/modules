process FASTQDL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed693e48ea3639b0a0f6a0ad191ad31f4e8fa15674a42c4053df6e06eff30cdc/data'
:         'community.wave.seqera.io/library/fastq-dl:4.0.1--0758300a069ba945' }"

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
