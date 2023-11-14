process BEDGOVCF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedgovcf:0.1.0--h9ee0642_0':
        'biocontainers/bedgovcf:0.1.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bed), path(config)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedgovcf: \$(echo \$(bedgovcf --version 2>&1) | sed 's/^bedgovcf version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedgovcf: \$(echo \$(bedgovcf --version 2>&1) | sed 's/^bedgovcf version //' )
    END_VERSIONS
    """
}
