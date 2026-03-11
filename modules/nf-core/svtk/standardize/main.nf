process SVTK_STANDARDIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtk:0.0.20190615--py37h73a75cf_2':
        'biocontainers/svtk:0.0.20190615--py37h73a75cf_2' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path (fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('svtk'), val('0.0.20190615'), topic: versions, emit: versions_svtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def contigs = fai ? "--contigs ${fai}" : ""

    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    svtk standardize \\
        $args \\
        ${contigs} \\
        ${input} \\
        ${prefix}.vcf.gz \\
        $meta.caller

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    echo | gzip > ${prefix}.vcf.gz

    """
}
