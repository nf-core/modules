process MEMOTE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/memote:0.17.0--pyhdfd78af_0' :
        'quay.io/biocontainers/memote:0.17.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(model)

    output:
    tuple val(meta), path("*.json.gz"), emit: json   , topic: json
    tuple val("${task.process}"), val('memote'), eval("memote --version | sed 's/memote, version //'"), topic: versions, emit: versions_memote

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\${PWD}
    export COBRA_SOLVER=glpk_exact

    memote run \\
        --filename ${prefix}.json.gz \\
        --ignore-git \\
        ${args} \\
        ${model}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\${PWD}
    export COBRA_SOLVER=glpk_exact
    
    echo "" | gzip > ${prefix}.json.gz
    """
}
