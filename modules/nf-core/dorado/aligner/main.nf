process DORADO_ALIGNER {
    tag "$meta.id"
    label 'process_high'

    // dorado is not on bioconda (ONTPL licence). Using
    // Docker Hub image directly. SHA tag pins to v1.4.0; a semver tag is tracked in
    // nanoporetech/dorado#1584. Same pattern as nf-core/parabricks modules.
    conda null
    container "docker.io/nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"

    input:
    tuple val(meta), path(bam)                          // unaligned BAM from dorado basecaller
    tuple val(meta2), path(reference), path(fai)        // reference FASTA (or .mmi index) and .fai

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.sequencing_summary.txt"), emit: summary, optional: true
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>&1 | head -1"), emit: versions_dorado, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
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
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def summary = args.contains('--emit-summary') ? "touch ${prefix}.sequencing_summary.txt" : ''
    """
    touch ${prefix}.bam
    ${summary}
    """
}
