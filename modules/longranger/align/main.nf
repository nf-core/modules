process LONGRANGER_ALIGN {
    // To use in cluster mode, some extra configurations is needed.
    // Visit tests/modules/longranger/align/nextflow.config for an example.

    tag "$meta.id"
    label 'process_medium'

    def version = '2.2.2-c1'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using longranger. Please use docker or singularity containers."
    }
    container "gitlab-registry.internal.sanger.ac.uk/tol-it/software/docker-images/longranger:${version}"

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
        longranger: \$(longranger align --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """
}
