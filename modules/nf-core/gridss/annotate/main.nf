process GRIDSS_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0'
        : 'biocontainers/gridss:2.13.2--h270b39a_0'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('gridss'), eval("GeneratePonBedpe --version 2>&1 | sed 's/-gridss//'"), topic: versions, emit: versions_gridss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.annotated"
    if ("${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    mkdir -p ${prefix}_work

    gridss_annotate_vcf_repeatmasker \\
        --output ${prefix}.vcf.gz \\
        --workingdir ${prefix}_work \\
        --threads ${task.cpus} \\
        ${args} \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.annotated"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
