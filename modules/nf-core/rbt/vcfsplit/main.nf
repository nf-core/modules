// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RBT_VCFSPLIT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e203258c4e52fcbe7f12392baf1bd184496aaeba52ceef17670924fa009245b1/data'
        : 'community.wave.seqera.io/library/rust-bio-tools:0.42.2--996dcdf3f4b075e4'}"

    input:
    tuple val(meta), path(vcf)
    val numchunks

    output:
    tuple val(meta), path("*.bcf"), emit: bcfchunks
    tuple val("${task.process}"), val('rbt'), eval("rbt --version | grep -oE '[0-9]+(\\.[0-9]+)+' | head -n 1"), topic: versions, emit: versions_rbt

    when:
    task.ext.when == null || task.ext.when

    script:
    // No args because tool does not accept args, only stdin/stdout
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunks = numchunks ? (numchunks - 1) : 15
    """
    rbt vcf-split \\
        ${vcf} \\
        ${prefix}.{0..${chunks}}.bcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunks = numchunks ? (numchunks - 1) : 15
    def bcf_files = (0..chunks).collect { items -> "${prefix}.${items}.bcf" }.join(' ')
    """
    touch ${bcf_files}
    """
}
