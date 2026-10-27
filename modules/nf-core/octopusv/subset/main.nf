process OCTOPUSV_SUBSET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/octopusv:1.0.0--pyhdfd78af_0':
        'quay.io/biocontainers/octopusv:1.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(svcf_in)
    tuple val(meta2), path(sample_file)
    tuple val(meta3), path(caller_file)

    output:
    tuple val(meta), path("*.svcf"), emit: svcf
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_file_arg = sample_file ? "--sample-file ${sample_file}" : ''
    def caller_file_arg = caller_file ? "--caller-file ${caller_file}" : ''
    def json_summary = args.contains('--json-summary') ? "> ${prefix}.json" : ""
    """
    octopusv subset \\
        --input-file ${svcf_in} \\
        --output-file ${prefix}.svcf \\
        ${sample_file_arg} \\
        ${caller_file_arg} \\
        ${json_summary} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_file_arg = sample_file ? "--sample-file ${sample_file}" : ''
    def caller_file_arg = caller_file ? "--caller-file ${caller_file}" : ''
    def touch_json = args.contains('--json-summary') ? "touch ${prefix}.json" : ""
    """
    echo ${sample_file_arg} ${caller_file_arg} $args

    touch ${prefix}.svcf
    ${touch_json}
    """
}
