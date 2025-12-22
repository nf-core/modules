process GENMOD_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmod:3.10.2--pyh7e72e81_0':
        'biocontainers/genmod:3.10.2--pyh7e72e81_0' }"

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
