

process DYSGU {
    tag "$meta.id"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dysgu:latest':
        'biocontainers/dysgu:latest' }"

    input:

    tuple val(meta), path(input_bam), path(input_bam_index)
    tuple val(meta2), path(reference_fasta), path(reference_fasta.fai)
    tuple val(meta3), path(temp_dir)

    output:

    tuple val(meta), path("*.vcf"), emit: vcf
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    dysgu run \\
        $reference \\
        $temp_dir \\
        $input_bam \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dysgu: \$(dysgu --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dysgu: \$(dysgu --version 2>&1)
    END_VERSIONS
    """
}
