process MIXCR_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    secret 'MI_LICENSE'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/milaboratory/mixcr/mixcr:4.7.0':
        'ghcr.io/milaboratory/mixcr/mixcr:4.7.0' }"

    input:
    tuple val(meta), path(reads)
    val preset
    val species

    output:
    tuple val(meta), path("*clones*.tsv"), emit: clones
    tuple val(meta), path("*.txt"),        emit: reports
    tuple val(meta), path("*.clns"),       emit: clns, optional: true
    tuple val(meta), path("*.vdjca"),      emit: vdjca, optional: true
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def java_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ''
    """
    mixcr $java_mem analyze \\
        $preset \\
        --species $species \\
        $args \\
        --threads $task.cpus \\
        ${reads} \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr -v 2>&1 | sed -n '1p' | sed -E 's/MiXCR v([0-9\\.]+).*/\1/' || true)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mixcr -v 
    touch ${prefix}.clones_TRA.tsv
    touch ${prefix}.clones_TRB.tsv
    touch ${prefix}.clns
    touch ${prefix}.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr -v 2>&1 | sed -n '1p' | sed -E 's/MiXCR v([0-9\\.]+).*/\1/' || true)
    END_VERSIONS
    """
}
