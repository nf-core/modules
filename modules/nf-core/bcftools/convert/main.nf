process BCFTOOLS_CONVERT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    path bed

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf_gz, optional: true
    tuple val(meta), path("*.vcf"), emit: vcf, optional: true
    tuple val(meta), path("*.bcf.gz"), emit: bcf_gz, optional: true
    tuple val(meta), path("*.bcf"), emit: bcf, optional: true
    tuple val(meta), path("*.hap.gz"), emit: hap, optional: true
    tuple val(meta), path("*.legend.gz"), emit: legend, optional: true
    tuple val(meta), path("*.samples"), emit: samples, optional: true
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def regions = bed ? "--regions-file ${bed}" : ""
    def reference = fasta ? "--fasta-ref ${fasta}" : ""
    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : args.contains("--haplegendsample") || args.contains("-h")
                        ? ""
                        : "vcf.gz"

    def output_cmd = args.contains("--haplegendsample") ? "" : "--output ${prefix}.${extension}"

    if ("${input}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    """
    bcftools convert \\
        ${args} \\
        ${regions} \\
        ${output_cmd} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : args.contains("--haplegendsample") || args.contains("-h")
                        ? "hap.gz"
                        : "vcf.gz"

    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""

    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    if ("${input}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    def hap = ""

    if (args.contains('--haplegendsample')) {
        def args_split = args.split(' ')
        hap = args_split.findIndexOf { arg -> arg == '--haplegendsample' }
        prefix = args_split[hap + 1]
    }
    """
    if [ -n "${hap}" ] ; then
        ${create_cmd} ${prefix}.hap.gz
        ${create_cmd} ${prefix}.legend.gz
        touch ${prefix}.samples
    else
        ${create_cmd} ${prefix}.${extension}
        ${create_index}
    fi
    """
}
