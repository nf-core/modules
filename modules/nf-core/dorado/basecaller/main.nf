// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DORADO_BASECALLER {
    tag "${meta.id}"
    // process_gpu only requests the accelerator; process_low/process_long set the
    // CPU, memory and time budget for the surrounding task.
    label 'process_low'
    label 'process_long'
    label 'process_gpu'

    // dorado is not on bioconda (ONTPL licence). Using
    // Docker Hub image directly. SHA tag pins to v2.1.1
    container "docker.io/nanoporetech/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(pod5)
    val model
    tuple val(meta2), path(models_dir)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
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
    // basecall on the GPU only when the pipeline requests one
    def device = task.accelerator ? 'cuda:all' : 'cpu'
    def models_arg = models_dir ? "--models-directory ${models_dir}" : "--models-directory ."
    """
    dorado \\
        basecaller \\
        ${args} \\
        --device ${device} \\
        ${models_arg} \\
        ${model} \\
        ${pod5} \\
        > ${prefix}.bam \\
        2>| >( tee ${prefix}.log >&2 )

    # --emit-summary writes a fixed-name `sequencing_summary.txt` into the working
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
    touch ${prefix}.log
    ${summary}
    """
}
