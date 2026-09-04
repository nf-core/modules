process EASYPQP_LIBRARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/easypqp:0.1.59--pyhdfd78af_0'
        : 'quay.io/biocontainers/easypqp:0.1.59--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(psmpkl), path(peakpkl)

    output:
    tuple val(meta), path("${prefix}.tsv")   , emit: tsv
    tuple val(meta), path("*_run_peaks.tsv"), emit: run_peaks, optional: true
    tuple val("${task.process}"), val('easypqp'), eval("easypqp --version 2>&1 | sed -nE 's/.*version ([0-9.]+).*/\\1/p'"), topic: versions, emit: versions_easypqp

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p \$MPLCONFIGDIR \$XDG_CACHE_HOME

    easypqp library \\
        --out ${prefix}.tsv \\
        ${args} \\
        ${psmpkl} ${peakpkl}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
