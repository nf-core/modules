process MINIMAC4_IMPUTE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/minimac4:4.1.6--hcb620b3_1'
        : 'quay.io/biocontainers/minimac4:4.1.6--hcb620b3_1'}"

    input:
    tuple val(meta), path(target_vcf), path(target_index), path(ref_msav), path(map), val(region)
    val write_sites

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: vcf
    tuple val(meta), path("${prefix}.sites.vcf.gz"), emit: sites
    tuple val("${task.process}"), val('minimac4'), eval("minimac4 --version |& sed '1!d ; s/minimac v//'"), emit: versions_minimac4, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args   ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-format bcf")    || args.contains("-O bcf")    ? "bcf"    :
                args.contains("--output-format sav")    || args.contains("-O sav")    ? "sav"    :
                args.contains("--output-format vcf.gz") || args.contains("-O vcf.gz") ? "vcf.gz" :
                args.contains("--output-format vcf")    || args.contains("-O vcf")    ? "vcf"    :
                args.contains("--output-format ubcf")   || args.contains("-O ubcf")   ? "ubcf"   :
                args.contains("--output-format usav")   || args.contains("-O usav")   ? "usav"   :
                "vcf.gz"
    def map_cmd      = map         ? "--map ${map}"                   : ""
    def region_cmd   = region      ? "--region ${region}"             : ""
    def sites_output = write_sites ? "--sites ${prefix}.sites.vcf.gz" : ""
    """
    minimac4 \\
        ${ref_msav} \\
        ${target_vcf} \\
        ${args} \\
        ${map_cmd} \\
        ${region_cmd} \\
        ${sites_output} \\
        --threads ${task.cpus} \\
        -o ${prefix}.${extension}
    """

    stub:
    def args  = task.ext.args   ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-format bcf")    || args.contains("-O bcf")    ? "bcf"    :
                args.contains("--output-format sav")    || args.contains("-O sav")    ? "sav"    :
                args.contains("--output-format vcf.gz") || args.contains("-O vcf.gz") ? "vcf.gz" :
                args.contains("--output-format vcf")    || args.contains("-O vcf")    ? "vcf"    :
                args.contains("--output-format ubcf")   || args.contains("-O ubcf")   ? "ubcf"   :
                args.contains("--output-format usav")   || args.contains("-O usav")   ? "usav"   :
                "vcf.gz"
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def sites_output = write_sites ? "echo '' | gzip > ${prefix}.sites.vcf.gz" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${sites_output}
    """
}
