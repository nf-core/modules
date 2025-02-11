process PYPGX_CREATEINPUTVCF {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    val(pgx_genes)
    val(assembly_version)



    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def assembly = "${assembly_version}" ?: "GRCh38"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pgx_genes = "--genes ${pgx_genes.join(' ')}" ?: ''

    """
    pypgx create-input-vcf \\
        ${args} \\
        ${pgx_genes} \\
        --assembly ${assembly} \\
        ${prefix}.vcf.gz \\
        ${fasta} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
