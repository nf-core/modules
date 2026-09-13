process CRABS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.0.7--pyhdfd78af_0' :
        'quay.io/biocontainers/crabs:1.0.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(crabsdb, stageAs: 'input/*', arity: '1..*')

    output:
    tuple val(meta), path("*.txt"), emit: crabsdb
    tuple val("${task.process}"), val('crabs'), eval("crabs --help 2>/dev/null | grep -oE 'v[0-9.]+' | cut -c2-"), emit: versions_crabs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = crabsdb.join(';')
    """
    crabs --merge \\
        --input '${input_files}' \\
        --output ${prefix}.txt \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
