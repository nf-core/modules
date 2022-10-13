process PLINK_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::plink=1.90b6.21" : null)
    def container_image = "/plink:1.90b6.21--h779adbc_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
}
