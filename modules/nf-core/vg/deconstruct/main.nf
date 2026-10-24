// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VG_DECONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc480ae940534804b86db6edfcbd5ed096d06443ba8062b7808d297cd37dfe6f/data'
:         'community.wave.seqera.io/library/vg:1.73.0--d74ddeb77d3e6db3' }"

    input:
    tuple val(meta), path(gfa)
    path(pb)
    path(gbwt)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('vg'), eval("vg 2>&1 | sed -n 's/.*version v\\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_vg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def snarls = pb ? "--snarls ${pb}" : ""
    def gbwt_arg = gbwt ? "--gbwt ${gbwt}" : ""
    """
    vg deconstruct \\
        --threads $task.cpus \\
        $args \\
        $snarls \\
        $gbwt_arg \\
        $gfa > ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
