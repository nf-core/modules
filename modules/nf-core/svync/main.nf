process SVYNC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svync:0.3.0--h9ee0642_0':
        'biocontainers/svync:0.3.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(config)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    tuple val("${task.process}"), val('svync'), eval("svync --version | sed 's/svync version //'"), emit: versions_svync, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def args3  = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    svync \\
        $args \\
        --config $config \\
        --input $vcf \\
        | bgzip --threads $task.cpus $args2 > ${prefix}.vcf.gz \\
        && tabix $args3 ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    echo | gzip -n > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
