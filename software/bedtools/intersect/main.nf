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

    //container "docker.pkg.github.com/nf-core/$MODULE"
    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    conda "${moduleDir}/environment.yml"

    input:
    path (input_file_1)
    path (input_file_2)
    val (intersectbed_args)

    output:
    stdout()

    script:
    """
    bedtools intersect -a ${input_file_1} -b ${input_file_2} ${intersectbed_args}
    """
}

