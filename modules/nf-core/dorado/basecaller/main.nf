process DORADO_BASECALLER {
    tag "$meta.id"
    // process_gpu only requests the accelerator; process_low/process_long set the
    // CPU, memory and time budget for the surrounding task.
    label 'process_low'
    label 'process_long'
    label 'process_gpu'

    // dorado is not on bioconda (ONTPL licence). Using
    // Docker Hub image directly. SHA tag pins to v1.4.0; a semver tag is tracked in
    // nanoporetech/dorado#1584. Same pattern as nf-core/parabricks modules.
    conda null
    container "docker.io/nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"

    input:
    tuple val(meta), path(pod5)                         // pod5 file or directory of pod5 files
    val(model)                                          // combined model string e.g. "sup,5mCG_5hmCG@latest", "hac@v5.0.0"
    tuple val(meta2), path(models_dir)                  // optional pre-downloaded models directory; pass [[],[]] to auto-download
    tuple val(meta3), path(reference), path(fai)        // optional reference FASTA for alignment; pass [[],[],[]] to skip

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.sequencing_summary.txt"), emit: summary, optional: true
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>&1 | head -1"), emit: versions_dorado, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def models_arg   = models_dir ? "--models-directory ${models_dir}" : "--models-directory ."
    def ref_arg      = reference  ? "--reference ${reference}"          : ""

    """
    dorado \\
        basecaller \\
        ${args} \\
        --device cuda:all \\
        ${models_arg} \\
        ${ref_arg} \\
        ${model} \\
        ${pod5} \\
        > ${prefix}.bam

    # --emit-summary writes a fixed-name `sequencing_summary.txt` into the working
    # directory; prefix it so outputs stay per-sample.
    if [ -f sequencing_summary.txt ]; then
        mv sequencing_summary.txt ${prefix}.sequencing_summary.txt
    fi
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def summary = args.contains('--emit-summary') ? "touch ${prefix}.sequencing_summary.txt" : ''
    """
    touch ${prefix}.bam
    ${summary}
    """
}
