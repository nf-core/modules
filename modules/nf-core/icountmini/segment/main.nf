process ICOUNTMINI_SEGMENT {
    tag "$gtf"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:3.0.1--pyh7cba7a3_0':
        'quay.io/biocontainers/icount-mini:3.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(gtf)
    path fai

    output:
    tuple val(meta), path("*_seg.gtf")       ,  emit: gtf
    tuple val(meta), path("*_regions.gtf.gz"),  emit: regions
    tuple val("${task.process}"), val('iCount-Mini'), eval("iCount-Mini -v"), emit: versions_icount_mini, topic: versions


    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${gtf.simpleName}"
    """
    iCount-Mini segment \\
        $args \\
        $gtf \\
        ${prefix}_seg.gtf \\
        $fai

    mv regions.gtf.gz ${prefix}_regions.gtf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.simpleName}"
    """
    touch ${prefix}_seg.gtf
    echo "" | gzip > ${prefix}_regions.gtf.gz
    """
}
