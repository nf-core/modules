process BCFTOOLS_REHEADER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(vcf), path(header), path(samples)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.{csi,tbi}"), emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fai_argument = fai ? "--fai ${fai}" : ""
    def header_argument = header ? "--header ${header}" : ""
    def samples_argument = samples ? "--samples ${samples}" : ""

    def args2 = task.ext.args2 ?: '--output-type z'
    def extension = args2.contains("--output-type b") || args2.contains("-Ob")
        ? "bcf.gz"
        : args2.contains("--output-type u") || args2.contains("-Ou")
            ? "bcf"
            : args2.contains("--output-type z") || args2.contains("-Oz")
                ? "vcf.gz"
                : args2.contains("--output-type v") || args2.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    """
    bcftools \\
        reheader \\
        ${fai_argument} \\
        ${header_argument} \\
        ${samples_argument} \\
        ${args} \\
        --threads ${task.cpus} \\
        ${vcf} \\
        | bcftools view \\
        ${args2} \\
        --output ${prefix}.${extension}
    """

    stub:
    def args2 = task.ext.args2 ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args2.contains("--output-type b") || args2.contains("-Ob")
        ? "bcf.gz"
        : args2.contains("--output-type u") || args2.contains("-Ou")
            ? "bcf"
            : args2.contains("--output-type z") || args2.contains("-Oz")
                ? "vcf.gz"
                : args2.contains("--output-type v") || args2.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index = args2.contains("--write-index=tbi") || args2.contains("-W=tbi")
        ? "tbi"
        : args2.contains("--write-index=csi") || args2.contains("-W=csi")
            ? "csi"
            : args2.contains("--write-index") || args2.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
