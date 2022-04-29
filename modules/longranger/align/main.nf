process LONGRANGER_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using longranger"
    }
    if ( workflow.containerEngine == 'singularity' || \
            workflow.containerEngine == 'docker' ) {
        exit 1, "Longranger can not be run in container environment"
    }

    input:
    tuple val(meta), path(reference)
    path(fastqs)

    output:
    tuple val(meta), path("${meta.id}/outs/possorted_bam.bam"), emit: bam
    tuple val(meta), path("${meta.id}/outs/possorted_bam.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}/outs/summary.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def jobmode = task.executor ?: "local"
    def sample = "${meta.id}"
    """
    longranger align --id=$sample --fastqs=$fastqs \
        --sample=$sample --reference=$reference \
        --jobmode=${jobmode} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    longranger: \$(echo \$(longranger mkref --version) | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """
}
