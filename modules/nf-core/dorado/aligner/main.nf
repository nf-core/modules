// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DORADO_ALIGNER {
    tag "${meta.id}"
    label 'process_high'

    // dorado is not on bioconda (ONTPL licence). Using
    // Docker Hub image directly. SHA tag pins to v2.1.1
    container "docker.io/nanoporetech/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(reference), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.sequencing_summary.txt"), emit: summary, optional: true
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>&1 | head -1"), emit: versions_dorado, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Dorado module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado \\
        aligner \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${bam} \\
        > ${prefix}.bam

    # --emit-summary writes a fixed-name `sequencing_summary.txt` into the work
    # directory; prefix it so outputs stay per-sample.
    if [ -f sequencing_summary.txt ]; then
        mv sequencing_summary.txt ${prefix}.sequencing_summary.txt
    fi
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Dorado module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def summary = args.contains('--emit-summary') ? "touch ${prefix}.sequencing_summary.txt" : ''
    """
    touch ${prefix}.bam
    ${summary}
    """
}
