// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VIENNARNA_RNACOFOLD {
    tag "$rnacofold_fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/1048b2d5571e2deab1e43e6629bbe4e6c7e9d8780b9dfbf0d3e0f44adb131021/data'
        : 'community.wave.seqera.io/library/viennarna:2.7.2--6e8b02540b675a6c'}"

    input:
    tuple val(meta), path(rnacofold_fasta)

    output:
    tuple val(meta), path("*.csv")                                                                   , emit: rnacofold_csv
    tuple val(meta), path("*.ps")                                                                    , emit: rnacofold_ps
    tuple val("${task.process}"), val('RNAcofold'), eval("RNAcofold --version 2>&1 | sed -n '1s/RNAcofold //p'"), emit: versions_rnacofold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNAcofold < ${rnacofold_fasta} \\
        --jobs=${task.cpus} \\
        --output-format="D" \\
        ${args} \\
        > ${rnacofold_fasta.baseName}.csv
    """

    stub:
    """
    touch ${rnacofold_fasta.baseName}.csv
    touch ${rnacofold_fasta.baseName}.ps
    """
}
