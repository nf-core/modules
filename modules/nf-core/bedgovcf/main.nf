process BEDGOVCF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedgovcf:0.1.1--h9ee0642_1':
        'biocontainers/bedgovcf:0.1.1--h9ee0642_1' }"

    input:
    tuple val(meta), path(bed), path(config)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val("bedgovcf"), eval("bedgovcf --version 2>&1 | sed 's/^bedgovcf version //'"), emit: versions_bedgovcf, topic: versions
    tuple val("${task.process}"), val("bgzip"), eval('bgzip --version | head -1 | sed "s/bgzip (htslib) //"')    , emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedgovcf \\
        $args \\
        --bed $bed \\
        --fai $fai \\
        --config $config \\
        | bgzip --stdout --threads $task.cpus $args2 > ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
