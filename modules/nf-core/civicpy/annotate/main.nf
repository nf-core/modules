process CIVICPY_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/civicpy:5.2.0--pyhdfd78af_0'
        : 'docker.io/griffithlab/civicpy:v5.2.0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val annotation_genome_version

    output:
    tuple val(meta), path("*.vcf"),                                                                                   emit: vcf
    tuple val("${task.process}"), val('civicpy'), eval("civicpy --version | sed 's/.*version //'"), topic: versions, emit: versions_civicpy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.civic"

    """
    export CIVICPY_CACHE_FILE=\$PWD/.civicpy

    civicpy annotate-vcf \\
        --input-vcf ${vcf} \\
        --output-vcf ${prefix}.vcf \\
        --reference ${annotation_genome_version} \\
        ${args}
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.civic"

    """
    touch ${prefix}.vcf
    """
}
