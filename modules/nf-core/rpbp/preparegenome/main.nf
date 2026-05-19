process RPBP_PREPAREGENOME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(fasta), path(gtf), path(config_yaml)

    output:
    tuple val(meta), path("rpbp_index"), path(config_yaml)                            , emit: index
    tuple val(meta), path("rpbp_index/*.annotated.bed.gz")                            , emit: transcript_bed
    tuple val(meta), path("rpbp_index/transcript-index/*.orfs-genomic.bed.gz")        , emit: orfs_genomic_bed
    tuple val(meta), path("rpbp_index/transcript-index/*.orfs-exons.bed.gz")          , emit: orfs_exons_bed
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mem_arg = task.memory ? "--mem ${task.memory.toGiga()}G" : ''
    """
    mkdir -p rpbp_index
    prepare-rpbp-genome \\
        ${config_yaml} \\
        --num-cpus ${task.cpus} \\
        ${mem_arg} \\
        ${args}
    """

    stub:
    """
    mkdir -p rpbp_index/transcript-index
    echo | gzip > rpbp_index/reference.annotated.bed.gz
    echo | gzip > rpbp_index/transcript-index/reference.orfs-genomic.bed.gz
    echo | gzip > rpbp_index/transcript-index/reference.orfs-exons.bed.gz
    """
}
