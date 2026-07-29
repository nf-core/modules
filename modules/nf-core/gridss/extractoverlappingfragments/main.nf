process GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3':
        'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(target_bed)

    output:
    tuple val(meta), path("*.subset.bam"),       emit: bam
    tuple val("${task.process}"), val('gridss'), eval("CallVariants --version 2>&1 | sed 's/-gridss//'"), topic: versions, emit: versions_gridss

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.subset"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    gridss_extract_overlapping_fragments \\
        -w '.' \\
        --targetbed  ${target_bed}  \\
        -o ${prefix}.bam \\
        $bam
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.subset"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    """
}
