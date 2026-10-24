process RIBOTRICER_PREPAREORFS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f634421ff566b498d7c9d065bdcce907e7cb62022f519417cabea6c3a921b41d/data'
        : 'community.wave.seqera.io/library/ribotricer_python_scipy:5721a090d55db3f6'}"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("*_candidate_orfs.tsv"), emit: candidate_orfs
    tuple val("${task.process}"), val('ribotricer'), eval("ribotricer --version 2>&1 | sed -n 's/^ribotricer, version //p'"), topic: versions, emit: versions_ribotricer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ribotricer prepare-orfs \\
        --gtf ${gtf} \\
        --fasta ${fasta} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_candidate_orfs.tsv
    """
}
