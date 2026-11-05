// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VIENNARNA_RNALFOLD {
    tag '$rna_fastq'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/1048b2d5571e2deab1e43e6629bbe4e6c7e9d8780b9dfbf0d3e0f44adb131021/data'
        : 'community.wave.seqera.io/library/viennarna:2.7.2--6e8b02540b675a6c'}"

    input:
    path fasta

    output:
    path "*.lfold"                                                                       , emit: rnalfold_txt
    tuple val("${task.process}"), val('RNALfold'), eval("RNALfold --version 2>&1 | sed -n '1s/RNALfold //p'"), emit: versions_rnalfold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNALfold \\
        ${args} \\
        --infile=$fasta \\
        --outfile=${fasta.baseName}.lfold
    """

    stub:
    """
    touch ${fasta.baseName}.lfold
    touch ${fasta}.ps
    """
}
