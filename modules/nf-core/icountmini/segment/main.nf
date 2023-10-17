process ICOUNTMINI_SEGMENT {
    tag "$gtf"
    label "process_single"

    conda "bioconda::icount-mini=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0' :
        'quay.io/biocontainers/icount-mini:2.0.3--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(gtf)
    path fai

    output:
    tuple val(meta), path("*_seg.gtf")         ,  emit: gtf
    tuple val(meta), path("*_regions.gtf.gz")  ,  emit: regions
    path "versions.yml"                        ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.simpleName}_seg"
    def regions_prefix = task.ext.regions_prefix ?: "${gtf.simpleName}"
    """
    iCount-Mini segment \\
        $gtf \\
        ${prefix}.gtf \\
        $fai

    mv regions.gtf.gz ${regions_prefix}_regions.gtf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}
