def MODULE = "bedtools_intersect"
params.publish_dir = MODULE
params.publish_results = "default"

process INTERSECT_BED {
    tag "$input_file_1-$input_file_2"

    publishDir "${params.out_dir}/${params.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (params.publish_results == "none") null
                    else filename }

    container "docker.pkg.github.com/nf-core/$MODULE"

    conda "${moduleDir}/environment.yml"

    input:
    path (input_file_1)
    path (input_file_2)
    val (intersectbed_args)

    output:
    path "${input_file_1.baseName}_i_${input_file_2.baseName}.bed", emit: intersect
    path "*.version.txt", emit: version

    script:
    def params_string = intersectbed_args.collect {
                            /-$it.key $it.value/
                        } join " "

    """
    bedtools intersect -a ${input_file_1} -b ${input_file_2} ${params_string} > ${input_file_1.baseName}_i_${input_file_2.baseName}.bed
    bedtools --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > bedtools.version.txt
    """
}
