// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MAFFT_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/355faf61033f2246b1f9f14e5f5a7e394fa5c08f32159986281686296d20affc/data'
        : 'community.wave.seqera.io/library/mafft:7.525--f4970860dd4f8697'}"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    tuple val("${task.process}"), val("mafft"), eval("mafft --version 2>&1 | sed 's/v//;s/ .*//'"), topic: versions, emit: versions_mafft

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mafft \\
        --thread ${task.cpus} \\
        --retree 0 \\
        ${args} \\
        --treeout ${fasta} > log.txt

    mv *.tree ${prefix}.dnd.tmp

    # remove all prefixes added by mafft which make the output incompatible with other tools
    awk '{gsub(/^[0-9]+_/, ""); print}' ${prefix}.dnd.tmp > ${prefix}.dnd
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dnd
    """
}
