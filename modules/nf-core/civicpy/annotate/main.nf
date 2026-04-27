process CIVICPY_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/civicpy:5.3.0--pyhdfd78af_0'
        : 'biocontainers/civicpy:5.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val annotation_genome_version
    path cache

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('civicpy'), eval("civicpy --version | sed 's/.*version //'"), topic: versions, emit: versions_civicpy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    export CIVICPY_CACHE_FILE=\$PWD/${cache}

    civicpy annotate-vcf \\
        --input-vcf ${vcf} \\
        --output-vcf ${prefix}.vcf \\
        --reference ${annotation_genome_version} \\
        ${args}

    bgzip ${prefix}.vcf
    """

    stub:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
