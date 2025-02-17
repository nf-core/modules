process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delly:1.3.3--h4d20210_0' :
        'biocontainers/delly:1.3.3--h4d20210_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(vcf), path(vcf_index), path(exclude_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.{bcf,vcf.gz}")  , emit: bcf
    tuple val(meta), path("*.{csi,tbi}")     , emit: csi
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"

    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""

    def bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    def vcf_output = suffix == "vcf" ? "| bgzip ${args2} --threads ${task.cpus} --stdout > ${prefix}.vcf.gz && tabix ${prefix}.vcf.gz" : ""

    def genotype = vcf ? "--vcffile ${vcf}" : ""

    """
    delly \\
        call \\
        ${args} \\
        ${bcf_output} \\
        --genome ${fasta} \\
        ${genotype} \\
        ${exclude} \\
        ${input} \\
        ${vcf_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"

    def bcf_output = suffix == "bcf" ? "touch ${prefix}.bcf && touch ${prefix}.bcf.csi" : ""
    def vcf_output = suffix == "vcf" ? "echo '' | gzip > ${prefix}.vcf.gz && touch ${prefix}.vcf.gz.tbi" : ""

    """
    ${bcf_output}
    ${vcf_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
