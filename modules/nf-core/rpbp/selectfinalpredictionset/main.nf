process RPBP_SELECTFINALPREDICTIONSET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta),  path(bayes_factors)
    tuple val(meta2), path(genome_fasta)

    output:
    tuple val(meta), path("${prefix}.predicted-orfs.bed.gz")     , emit: predicted
    tuple val(meta), path("${prefix}.predicted-orfs.dna.fa")     , emit: dna_fasta
    tuple val(meta), path("${prefix}.predicted-orfs.protein.fa") , emit: protein_fasta
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
        ${prefix}.predicted-orfs.bed.gz \\
        ${prefix}.predicted-orfs.dna.fa \\
        ${prefix}.predicted-orfs.protein.fa \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.predicted-orfs.bed.gz
    touch ${prefix}.predicted-orfs.dna.fa
    touch ${prefix}.predicted-orfs.protein.fa
    """
}
