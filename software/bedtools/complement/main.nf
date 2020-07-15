def MODULE = "bedtools_complement"
params.publish_dir = MODULE
params.publish_results = "default"

process BEDTOOLS_COMPLEMENT {
    tag {input_file}

    publishDir "${params.out_dir}/${params.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (params.publish_results == "none") null
                    else filename }

    //container "docker.pkg.github.com/nf-core/$MODULE"
    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    conda "${moduleDir}/environment.yml"

    input:
    path (input_file)
    path (fasta_sizes)
    val (bedtools_complement_args)

    output:
    path "${input_file}.bed", emit: complement
    path "*.version.txt", emit: version

    script:
    """
    bedtools complement -i ${input_file} -g ${fasta_sizes} ${bedtools_complement_args} > ${input_file}.bed
    bedtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > bedtools.version.txt
    """
}
