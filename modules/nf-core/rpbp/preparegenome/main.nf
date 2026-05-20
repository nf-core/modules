process RPBP_PREPAREGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("${prefix}")                                                       , emit: index
    tuple val(meta), path("${prefix}/${name}.annotated.bed.gz")                              , emit: transcript_bed
    tuple val(meta), path("${prefix}/transcript-index/${name}.orfs-genomic.annotated.bed.gz"), emit: orfs_genomic_bed
    tuple val(meta), path("${prefix}/transcript-index/${name}.orfs-exons.annotated.bed.gz")  , emit: orfs_exons_bed
    path "versions.yml"                                                                      , emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    name   = meta.id ?: 'reference'
    prefix = task.ext.prefix ?: 'rpbp_index'
    template 'prepare_rpbp_genome.py'

    stub:
    name   = meta.id ?: 'reference'
    prefix = task.ext.prefix ?: 'rpbp_index'
    """
    mkdir -p ${prefix}/transcript-index
    echo "" | gzip > ${prefix}/${name}.annotated.bed.gz
    echo "" | gzip > ${prefix}/transcript-index/${name}.orfs-genomic.annotated.bed.gz
    echo "" | gzip > ${prefix}/transcript-index/${name}.orfs-exons.annotated.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
        rpbp: \$(python -c "import rpbp; print(rpbp.__version__)")
    END_VERSIONS
    """
}
