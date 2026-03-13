process CNVPYTOR_VIEW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.2.1--pyhdfd78af_0':
        'biocontainers/cnvpytor:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pytor_files)
    val bin_sizes
    val output_format

    output:
    tuple val(meta), path("*.vcf"), emit: vcf      , optional: true
    tuple val(meta), path("*.tsv"), emit: tsv      , optional: true
    tuple val(meta), path("*.xls"), emit: xls      , optional: true
    tuple val("${task.process}"), val('cnvpytor'), eval("cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p'"), emit: versions_tool1, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_suffix = output_format ?: 'vcf'
    def bins   = bin_sizes ?: '1000'
    def input  = pytor_files.join(" ")
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    python3 <<CODE
    import cnvpytor,os
    binsizes = "${bins}".split(" ")
    for binsize in binsizes:
        file_list = "${input}".split(" ")
        app = cnvpytor.Viewer(file_list, params={} )
        outputfile = "{}_{}.{}".format("${prefix}",binsize.strip(),"${output_suffix}")
        app.print_filename = outputfile
        app.bin_size = int(binsize)
        app.print_calls_file()
    CODE
    """

    stub:
    def output_suffix = output_format ?: 'vcf'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${output_suffix}
    """
}
