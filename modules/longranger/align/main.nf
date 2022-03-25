
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
 
    label 'mem_high'

    input:
    tuple val(meta), val(sample)
    tuple val(meta), path(fastqs)
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path("${sample}/outs/possorted_bam.bam"), emit: bam
    tuple val(meta), path("${sample}/outs/possorted_bam.bam.bai"), emit: bai
    tuple val(meta), path("${sample}/outs/summary.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${moduleDir}/bin/run_longranger.py $sample $fastqs $reference ${moduleDir}/override.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    longranger align --version
    END_VERSIONS
    """
}
