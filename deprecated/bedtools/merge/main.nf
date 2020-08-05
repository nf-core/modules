def MODULE = "bedtools_merge"
params.publish_dir = MODULE
params.publish_results = "default"

process BEDTOOLS_MERGE {
    tag { input_file }

    publishDir "${params.out_dir}/${params.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (params.publish_results == "none") null
                    else filename }

    container "docker.pkg.github.com/nf-core/$MODULE"

    conda "${moduleDir}/environment.yml"

    input:
    path (input_file)
    val (bedtools_merge_args)

    output:
    path "${input_file}.bed", emit: merge
    path "*.version.txt", emit: version

    script:
    """
    bedtools merge -i ${input_file} ${bedtools_merge_args} > ${input_file}.bed
    bedtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > bedtools.version.txt
    """
}
