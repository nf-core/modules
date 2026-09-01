process OPT_FLIP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64550ef193f98ea294d70e087a80e37bbf0c5c4a5920f0a22414ed8a11c32caa/data' :
        'community.wave.seqera.io/library/opt:0.0.1--f0b1e63f50e38ab1'}"

    input:
    tuple val(meta), path(probes_fasta)
    tuple val(meta2), path(ref_annot_gff), path(ref_annot_fa)

    output:
    tuple val(meta), path("${prefix}/fwd_oriented.fa"), emit: fwd_oriented_fa
    tuple val("${task.process}"), val('opt'), eval("opt --version"), topic: versions, emit: versions_opt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    opt \\
        -o ${prefix} \\
        -p ${task.cpus} \\
        flip \\
        -i ${probes_fasta} \\
        -a ${ref_annot_gff} \\
        -f ${ref_annot_fa} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/fwd_oriented.fa"
    """
}
