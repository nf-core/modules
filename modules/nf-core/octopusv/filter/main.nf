process OCTOPUSV_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/octopusv:1.0.0--pyhdfd78af_0':
        'quay.io/biocontainers/octopusv:1.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(svcf_in)

    output:
    tuple val(meta), path("*.svcf"), emit: svcf
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def json_summary = args.contains('--json-summary') ? "> ${prefix}.json" : ""
    """
    octopusv filter \\
        -i ${svcf_in} \\
        -o ${prefix}.svcf \\
        ${json_summary} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_json = args.contains('--json-summary') ? "touch ${prefix}.json" : ""
    """
    echo $args

    touch ${prefix}.svcf
    ${touch_json}
    """
}
