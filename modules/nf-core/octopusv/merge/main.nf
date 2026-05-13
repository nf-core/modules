process OCTOPUSV_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/octopusv:0.3.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/octopusv:0.3.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(svcfs)

    output:
    tuple val(meta), path("*.svcf"), emit: svcf
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = (task.ext.args ?: '').trim()
    def prefix = task.ext.prefix ?: "${meta.id}"

    def strategy_flags = ['--union', '--intersect', '--specific', '--min-support', '--exact-support', '--max-support', '--expression']

    def explicit_flags = strategy_flags.findAll { flag -> args.tokenize().contains(flag) }
    if (explicit_flags.size() > 1) {
        error("OCTOPUSV_MERGE: multiple merge strategies specified: ${explicit_flags.join(', ')}")
    }

    if (explicit_flags.isEmpty()) {
        args = "--union ${args}".trim()
    }

    """
    octopusv merge -i ${svcfs.join(' ')} \\
        -o ${prefix}.svcf \\
        ${args}
    """

    stub:
    def args = (task.ext.args ?: '').trim()
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.svcf
    """
}
