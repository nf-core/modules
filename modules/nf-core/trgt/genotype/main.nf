process TRGT_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trgt:1.5.0--h9ee0642_0':
        'biocontainers/trgt:1.5.0--h9ee0642_0' }"

    input:
    tuple val(meta) , path(bam), path(bai), val(karyotype)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(repeats)

    output:
    tuple val(meta), path("*.vcf.gz")      , emit: vcf
    tuple val(meta), path("*.spanning.bam"), emit: bam     , optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def karyo = karyotype ? "--karyotype ${karyotype}" : ""
    """
    trgt genotype \\
        $args \\
        --genome ${fasta} \\
        --reads ${bam} \\
        --repeats ${repeats} \\
        ${karyo} \\
        --threads ${task.cpus} \\
        --output-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$(trgt --version |& sed '1!d ; s/trgt //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spanning.bam
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$(trgt --version |& sed '1!d ; s/trgt //')
    END_VERSIONS
    """
}
