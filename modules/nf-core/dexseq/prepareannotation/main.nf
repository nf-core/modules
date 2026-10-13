process DEXSEQ_PREPAREANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f3e6fbbab7a2e4b6e1d57be3e4aa40b18810e3a72b7d6a8c0fcbb2bf0cc7c63/data'
        : 'community.wave.seqera.io/library/bioconductor-dexseq_htseq_python:dfa2c806cf6767c6'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val("${task.process}"), val('bioconductor-dexseq'), eval('Rscript -e "cat(as.character(packageVersion(\'DEXSeq\')))"'), topic: versions, emit: versions_dexseq
    tuple val("${task.process}"), val('htseq'), eval('python -c "import HTSeq; print(HTSeq.__version__)"'), topic: versions, emit: versions_htseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The script ships with the DEXSeq R package and is not on the PATH.
    // The script iterates over Python sets, so a fixed hash seed is needed to keep the
    // order of the aggregate gene IDs and of the transcripts attribute reproducible.
    // The order under a given seed depends on the Python version, which is pinned in environment.yml.
    """
    dexseq_scripts=\$(Rscript -e 'cat(system.file("python_scripts", package = "DEXSeq"))')

    PYTHONHASHSEED=0 python \${dexseq_scripts}/dexseq_prepare_annotation.py \\
        ${args} \\
        ${gtf} \\
        ${prefix}.gff
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff
    """
}
