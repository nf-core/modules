process CIVICPY_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7c82422a425f51e197c0f88b2927e5fdcf61afb752dec396a56511c21605931/data'
        : 'community.wave.seqera.io/library/civicpy_htslib:a4e0f1666f9ba596' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val annotation_genome_version
    path cache

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('civicpy'), eval("civicpy --version | sed 's/.*version //'"), topic: versions, emit: versions_civicpy
    tuple val("${task.process}"), val('htslib'),  eval("bgzip --version 2>&1 | head -1 | sed 's/bgzip (htslib) //'"),  topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    def cache_file = cache            ? "\$PWD/${cache}" : "\$PWD/${prefix}.civicpy_cache.pkl"
    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    export CIVICPY_CACHE_FILE=${cache_file}

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
