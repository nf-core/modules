process VCFLIB_VCFFIXUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/169e4e28f26469eb05baf60eab777bccadd747ac75038c6bb22149cd40c2ff38/data':
        'community.wave.seqera.io/library/bcftools_vcflib:0b47030679d1eff1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.{vcf,bcf}{,.gz}"), emit: vcf
    tuple val(meta), path("*.{csi,tbi}")      , emit: index, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('vcflib'), val("1.0.14"), topic: versions, emit: versions_vcflib
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools


    script:
    def args    = task.ext.args   ?: ''
    def args2   = task.ext.args2  ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}.fixed"

    def extension = args2.contains("--output-type b") || args2.contains("-Ob")
        ? "bcf.gz"
        : args2.contains("--output-type u") || args2.contains("-Ou")
            ? "bcf"
            : args2.contains("--output-type z") || args2.contains("-Oz")
                ? "vcf.gz"
                : args2.contains("--output-type v") || args2.contains("-Ov")
                    ? "vcf"
                    : "vcf"

    if ( "${vcf}" == "${prefix}.${extension}" ) {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    vcffixup \\
        ${args} \\
        ${vcf} \\
        | bcftools view ${args2} -o ${prefix}.${extension}
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}.fixed"
    def args2   = task.ext.args2  ?: ''
    def extension = args2.contains("--output-type b") || args2.contains("-Ob")
        ? "bcf.gz"
        : args2.contains("--output-type u") || args2.contains("-Ou")
            ? "bcf"
            : args2.contains("--output-type z") || args2.contains("-Oz")
                ? "vcf.gz"
                : args2.contains("--output-type v") || args2.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def stub_index = args2.contains("--write-index=tbi") || args2.contains("-W=tbi")
        ? "tbi"
        : args2.contains("--write-index=csi") || args2.contains("-W=csi")
            ? "csi"
            : args2.contains("--write-index") || args2.contains("-W")
                ? "csi"
                : ""

    if ( "${vcf}" == "${prefix}.${extension}" ) {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && stub_index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${stub_index}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
