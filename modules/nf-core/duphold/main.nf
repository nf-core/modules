process DUPHOLD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/duphold:0.2.1--h516909a_1':
        'biocontainers/duphold:0.2.1--h516909a_1' }"

    input:
    tuple val(meta), path(alignment_file), path(alignment_index), path(sv_variants), path(snp_variants), path(snp_variants_index)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def snp_annotation = snp_variants ? "--snp ${snp_variants}" : ""

    """
    duphold \\
        ${args} \\
        --threads ${task.cpus} \\
        --output ${prefix}.vcf.gz \\
        --vcf ${sv_variants} \\
        --bam ${alignment_file} \\
        --fasta ${fasta} \\
        ${snp_annotation}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duphold: \$(duphold -h | head -n 1 | sed -e "s/^version: //")
    END_VERSIONS
    """
}
