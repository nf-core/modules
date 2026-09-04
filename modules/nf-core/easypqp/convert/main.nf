process EASYPQP_CONVERT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/easypqp:0.1.59--pyhdfd78af_0'
        : 'quay.io/biocontainers/easypqp:0.1.59--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(id_file), path(spectra)
    path unimod

    output:
    tuple val(meta), path("*.psmpkl") , emit: psmpkl
    tuple val(meta), path("*.peakpkl"), emit: peakpkl
    tuple val("${task.process}"), val('easypqp'), eval("easypqp --version 2>&1 | sed -nE 's/.*version ([0-9.]+).*/\\1/p'"), topic: versions, emit: versions_easypqp

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def unimod_arg = unimod ? "--unimod ${unimod}" : ''
    // Output basenames are derived from the spectra file by easypqp and must not be
    // changed: easypqp library pairs peak and PSM pickles to runs by this basename.
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p \$MPLCONFIGDIR \$XDG_CACHE_HOME

    easypqp convert \\
        --pepxml ${id_file} \\
        --spectra ${spectra} \\
        ${unimod_arg} \\
        ${args}
    """

    stub:
    def prefix = spectra.baseName
    """
    touch ${prefix}.psmpkl
    touch ${prefix}.peakpkl
    """
}
