process BACPHLIP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f45b980d08c53e74f3b8687b03aa8569401c6079a444c0efe5e4096e1367a751/data'
        : 'community.wave.seqera.io/library/bacphlip_hmmer_numpy:23d5784a0d17b398' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bacphlip")         , emit: bacphlip_results
    tuple val(meta), path("*.hmmsearch.tsv")    , emit: hmmsearch_results
    tuple val("${task.process}"), val('bacphlip'), val('0.9.6'), emit: versions_bacphlip, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bacphlip \\
        -i $fasta \\
        $args
    """

    stub:
    """
    touch ${fasta}.bacphlip
    touch ${fasta}.hmmsearch.tsv
    """
}
