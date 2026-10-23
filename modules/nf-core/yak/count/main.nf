// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process YAK_COUNT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc0349f7b2bdd3d35e237ea545cb2523a91791cfa1246f4ac18c71aed303f2ab/data'
:         'community.wave.seqera.io/library/yak:0.1--3a5526a015428c21' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.yak"), emit: yak
    tuple val("${task.process}"), val('yak'), eval("yak version"), emit: versions_yak, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input_command = meta.single_end ? "${fastq}" : "<(zcat ${fastq}) <(zcat ${fastq})"
    """
    yak \\
        count \\
        $args \\
        -t${task.cpus} \\
        -o ${prefix}.yak \\
        $input_command
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yak
    """
}
