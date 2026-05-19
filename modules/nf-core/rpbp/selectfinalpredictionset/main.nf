process RPBP_SELECTFINALPREDICTIONSET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(bayes_factors)
    path  genome_fasta

    output:
    tuple val(meta), path("${prefix}.predicted-orfs.filtered.bed.gz")     , emit: predicted
    tuple val(meta), path("${prefix}.predicted-orfs.filtered.dna.fa")     , emit: dna_fasta
    tuple val(meta), path("${prefix}.predicted-orfs.filtered.protein.fa") , emit: protein_fasta
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    select-final-prediction-set \\
        ${bayes_factors} \\
        ${genome_fasta} \\
        ${prefix}.predicted-orfs.filtered.bed.gz \\
        ${prefix}.predicted-orfs.filtered.dna.fa \\
        ${prefix}.predicted-orfs.filtered.protein.fa \\
        --select-longest-by-stop \\
        --select-best-overlapping \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.predicted-orfs.filtered.bed.gz
    touch ${prefix}.predicted-orfs.filtered.dna.fa
    touch ${prefix}.predicted-orfs.filtered.protein.fa
    """
}
