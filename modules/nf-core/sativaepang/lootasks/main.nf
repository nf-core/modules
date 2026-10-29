process SATIVAEPANG_LOOTASKS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.10.0--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.10.0--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(refjson)

    output:
    tuple val(meta), path("*.l1o_tasks"), emit: taskdir
    tuple val("${task.process}"), val('sativaepang'), eval("sativa-epang --version | cut -d' ' -f2"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sativa-epang \\
        -stage loo-tasks \\
        -r ${refjson} \\
        -n ${prefix} \\
        -o . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.l1o_tasks/fold_000
    touch ${prefix}.l1o_tasks/manifest.json
    touch ${prefix}.l1o_tasks/fold_000/ref.nwk
    touch ${prefix}.l1o_tasks/fold_000/ref.fasta
    touch ${prefix}.l1o_tasks/fold_000/query.fasta
    """
}
