process VEMBRANE_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vembrane:2.4.0--pyhdfd78af_0'
        : 'biocontainers/vembrane:2.4.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf)
    val expression

    output:
    tuple val(meta), path("*.{vcf,bcf,bcf.gz}"), emit: vcf
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-fmt vcf") || args.contains("-Ovcf") ? "vcf" :
                    args.contains("--output-fmt bcf") || args.contains("-Obcf") ? "bcf" :
                    args.contains("--output-fmt uncompressed-bcf") || args.contains("-Ouncompressed-bcf") ? "bcf.gz" :
                    "vcf"

    if ("${vcf}" == "${prefix}.${extension}") {
        error("Input and output names are the same, use \"task.ext.prefix\" in module configuration to disambiguate!")
    }

    """
    vembrane filter \\
        ${args} \\
        ${expression} \\
        -o ${prefix}.${extension} \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vembrane: \$(vembrane --version | sed '1!d;s/.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-fmt vcf") || args.contains("-Ovcf") ? "vcf" :
                    args.contains("--output-fmt bcf") || args.contains("-Obcf") ? "bcf" :
                    args.contains("--output-fmt uncompressed-bcf") || args.contains("-Ouncompressed-bcf") ? "bcf.gz" :
                    "vcf"

    if ("${vcf}" == "${prefix}.${extension}") {
        error("Input and output names are the same, use \"task.ext.prefix\" in module configuration to disambiguate!")
    }

    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vembrane: \$(vembrane --version | sed '1!d;s/.* //')
    END_VERSIONS
    """
}
