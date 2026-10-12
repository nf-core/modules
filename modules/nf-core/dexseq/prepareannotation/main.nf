process DEXSEQ_PREPAREANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/0718f75da083443fe7acb310d3f8ef8216fb4561395226a48fff870d17323f4e/data'
        : 'community.wave.seqera.io/library/bioconductor-dexseq_htseq:4266bf085c9883ac'}"

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
    // A fixed hash seed keeps the order of the transcripts attribute reproducible.
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
