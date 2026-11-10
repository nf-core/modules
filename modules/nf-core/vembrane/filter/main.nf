process VEMBRANE_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6a/6a1bc5a7c85eebca118ae436bda77e005cc6c1ca0e60687185fb1d2f195e5eaf/data'
        : 'community.wave.seqera.io/library/vembrane:2.6.2--0bc429af0e13a402'}"

    input:
    tuple val(meta), path(vcf)
    val expression

    output:
    tuple val(meta), path("*.{vcf,bcf,bcf.gz}"), emit: vcf
    tuple val("${task.process}"), val('vembrane'), eval("vembrane --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_vembrane

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
    """
}
