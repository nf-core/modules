def MODULE = "bedtools_genomecov"
params.publish_dir = MODULE
params.publish_results = "default"

process BEDTOOLS_GENOMECOV {
    tag {bam}

    publishDir "${params.out_dir}/${params.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (params.publish_results == "none") null
                    else filename }

    //container "docker.pkg.github.com/nf-core/$MODULE"
    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    conda "${moduleDir}/environment.yml"

    input:
    path (bam)
    path (chrom_sizes)
    val (bedtools_genomecov_args)

    output:
    path "${bam}.bed", emit: coverage
    path "*.version.txt", emit: version

    script:
    """
    bedtools genomecov -ibam ${bam} -g ${chrom_sizes} ${bedtools_genomecov_args} > ${bam}.bed
    bedtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > bedtools.version.txt
    """
}
