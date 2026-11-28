// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METAMDBG_ASM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8cc5344c1450ae8914ecb76febc92ec2dbade6a79e301bd1c65e0133e5b5631f/data'
        : 'community.wave.seqera.io/library/metamdbg:1.4--553e8e26a319ac0d'}"

    input:
    tuple val(meta), path(reads, arity: '1..*')
    val(input_type)

    output:
    tuple val(meta), path("*.contigs.fasta.gz"), emit: contigs
    tuple val(meta), path("*.metaMDBG.log")    , emit: log
    tuple val("${task.process}"), val('metamdbg'), eval('metaMDBG | sed -n "s/.*Version: //p"'), emit: versions_metamdbg, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if(!(input_type in ["hifi", "ont"])) {
        error("ERROR: input_type must be one of either 'hifi' or 'ont'.")
    }
    """
    metaMDBG asm \\
        --threads ${task.cpus} \\
        --out-dir . \\
        ${args} \\
        --in-${input_type} ${reads}

    rm -r tmp/

    mv contigs.fasta.gz ${prefix}.contigs.fasta.gz
    mv metaMDBG.log ${prefix}.metaMDBG.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}.metaMDBG.log
    echo "" | gzip > ${prefix}.contigs.fasta.gz
    """
}
