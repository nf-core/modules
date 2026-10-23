// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VG_CONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc480ae940534804b86db6edfcbd5ed096d06443ba8062b7808d297cd37dfe6f/data'
:         'community.wave.seqera.io/library/vg:1.73.0--d74ddeb77d3e6db3' }"

    input:
    tuple val(meta), path(input), path(tbis), path(insertions_fasta)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.vg") , emit: graph
    tuple val("${task.process}"), val('vg'), eval("vg 2>&1 | sed -n 's/.*version v\\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_vg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def mode =  input instanceof ArrayList || input.toString().endsWith(".vcf.gz") ? 'vcf' : 'msa'

    input_files = mode == 'vcf' ? input.collect { vcf_file -> "--vcf ${vcf_file}" }.join(" ") : "--msa ${input}"
    reference = mode == 'vcf' ? "--reference ${fasta}" : ""

    insertions = insertions_fasta ? "--insertions ${insertions_fasta}" : ""

    """
    vg construct \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input_files} \\
        ${insertions} \\
        > ${prefix}.vg
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vg
    """
}
