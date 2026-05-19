process RPBP_EXTRACTMETAGENEPROFILES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  transcript_bed

    output:
    tuple val(meta), path("${prefix}.metagene-profile.csv.gz"), emit: metagene
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract-metagene-profiles \\
        ${bam} \\
        ${transcript_bed} \\
        ${prefix}.metagene-profile.csv.gz \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.metagene-profile.csv.gz
    """
}
