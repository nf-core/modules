process COUNTREADS {
    tag "$meta.id"
    label 'process_medium'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'drop':
        'drop' }"

    input:
        tuple val(meta), path(Rds)
        each  path(bam)

    output:
        tuple path("*.Rds")
        path "versions.yml"           , emit: versions


    when:
        task.ext.when == null || task.ext.when


    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /abberant-expression-pipeline/Counting/countReads.R $bam $Rds $csv $params

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aberrantexpression: \$(echo 1.3.0)
    END_VERSIONS
    """
}
