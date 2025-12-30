process GENMOD_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/31b331bee43c7ff070bdde5460a4102ba31c3bfb0ee0d70197001ff011036555/data' :
        'community.wave.seqera.io/library/genmod_python:31b2fba4d3b7ba6f' }"

    input:
    tuple val(meta), path(input_vcf)

    output:
    tuple val(meta), path("*_annotate.vcf"), emit: vcf
    tuple val("${task.process}"), val('genmod'), eval("genmod --version | sed 's/^.*genmod version: //'"), emit: versions_genmod, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genmod \\
        annotate \\
        $args \\
        --outfile ${prefix}_annotate.vcf \\
        $input_vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotate.vcf
    """
}
