process PREPROCESSGENEANNOTATION {
    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'drop':
        'drop' }"

    input:
        tuple val(meta), path(gtf)

    output:
        tuple path("txdb.db"), path("*.tsv"), path("count_ranges.Rds") , emit: count_ranges
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /Scripts/Pipeline/preprocessGeneAnnotation.R $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aberrantexpression: \$(echo 1.3.0)
    END_VERSIONS
    """
}
