process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_gpu'

    // No bioconda package for dorado (ONTPL licence, not on bioconda/conda-forge).
    // Docker image from nanoporetech/dorado Docker Hub — SHA tag maps to v1.4.0 (released 2026-02-19).
    // Track 2 (TODO): submit dorado to bioconda → auto-generate depot.galaxyproject.org Singularity image
    //                 and replace container spec with standard nf-core biocontainers format.
    conda null
    // Track 1 (lab use): local SIF — dorado not on bioconda, no depot.galaxyproject.org image available.
    // Track 2 (nf-core upstream): replace with biocontainers URI once dorado is on bioconda.
    //   Docker equivalent: nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e (= v1.4.0)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/data1/greenbab/software/images/dorado_compiled_v1.4.0.sif' :
        'nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e' }"

    input:
    tuple val(meta), path(pod5)                         // pod5 file or directory of pod5 files
    val(model)                                          // combined model string e.g. "sup,5mCG_5hmCG@latest", "hac@v5.0.0"
    tuple val(meta2), path(models_dir)                  // optional pre-downloaded models directory; pass [[],[]] to auto-download
    tuple val(meta3), path(reference), path(fai)        // optional reference FASTA for alignment; pass [[],[],[]] to skip

    output:
    tuple val(meta), path("*.bam")      , emit: bam
    tuple val(meta), path("*_summary.tsv"), emit: summary , optional: true
    tuple val(meta), path("*.log")      , emit: log       , optional: true
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>&1 | head -1 | sed 's/^//'"), emit: versions_dorado, topic: versions

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
        --device ${task.ext.device ?: 'cuda:all'} \\
        ${models_arg} \\
        ${ref_arg} \\
        ${model} \\
        ${pod5} \\
        > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}_summary.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: 1.4.0
    END_VERSIONS
    """
}
