process RPBP_EXTRACTORFPROFILES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(lengths_offsets)
    tuple val(meta2), path(orfs_genomic_bed)
    tuple val(meta3), path(exons_bed)

    output:
    tuple val(meta), path("${prefix}.profiles.mtx.gz"), emit: profiles
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    LENGTHS=\$(tail -n +2 ${lengths_offsets} | cut -f1 | tr '\\n' ' ')
    OFFSETS=\$(tail -n +2 ${lengths_offsets} | cut -f2 | tr '\\n' ' ')

    extract-orf-profiles \\
        ${bam} \\
        ${orfs_genomic_bed} \\
        ${exons_bed} \\
        ${prefix}.profiles.mtx.gz \\
        --lengths \$LENGTHS \\
        --offsets \$OFFSETS \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.profiles.mtx.gz
    """
}
