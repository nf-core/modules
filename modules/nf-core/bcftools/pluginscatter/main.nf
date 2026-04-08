process BCFTOOLS_PLUGINSCATTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val sites_per_chunk
    val scatter
    path scatter_file
    path regions
    path targets

    output:
    tuple val(meta), path("*{vcf,vcf.gz,bcf,bcf.gz}"), emit: scatter
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def mandatory_arg = sites_per_chunk ? "--nsites-per-chunk ${sites_per_chunk}" : scatter ? "--scatter ${scatter}" : "--scatter-file ${scatter_file}"
    def regions_arg = regions ? "--regions-file ${regions}" : ""
    def targets_arg = targets ? "--targets-file ${targets}" : ""
    """
    bcftools plugin scatter \\
        ${vcf} \\
        ${mandatory_arg} \\
        ${regions_arg} \\
        ${targets_arg} \\
        --output ${prefix} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    mv ${prefix}/* .
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
                    : "vcf"
    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index_1 = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}0.${extension}.${index}" : ""
    def create_index_2 = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}1.${extension}.${index}" : ""
    def create_index_3 = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}2.${extension}.${index}" : ""

    """
    ${create_cmd} ${prefix}0.${extension}
    ${create_cmd} ${prefix}1.${extension}
    ${create_cmd} ${prefix}2.${extension}

    ${create_index_1}
    ${create_index_2}
    ${create_index_3}
    """
}
