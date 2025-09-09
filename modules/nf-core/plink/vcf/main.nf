process PLINK_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    tuple val(meta), path("*.bim"), emit: bim, optional: true
    tuple val(meta), path("*.fam"), emit: fam, optional: true

    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink \\
        --vcf ${vcf} \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    touch ${prefix}.bim
    touch ${prefix}.fam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
