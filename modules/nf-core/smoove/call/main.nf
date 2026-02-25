process SMOOVE_CALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoove:0.2.8--h9ee0642_1' :
        'biocontainers/smoove:0.2.8--h9ee0642_1' }"

    input:
    tuple val(meta), path(input), path(index), path(exclude_beds)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('smoove'), eval("smoove -v |& sed -n 's/smoove version: *//p'"), emit: versions_smoove, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def exclude = exclude_beds ? "--exclude ${exclude_beds}" : ""
    """
    smoove call \\
        ${args} \\
        --outdir . \\
        --name ${prefix} \\
        --fasta ${fasta} \\
        ${exclude} \\
        --processes ${task.cpus} \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
