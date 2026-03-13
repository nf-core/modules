process CNVPYTOR_VIEW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bb/bbb6343edff4191cb1f445b2aac028d1f805ed5a7d50799513c82531bcfdede5/data':
        'community.wave.seqera.io/library/cnvpytor_make:a8fdcebe82041114' }"

    input:
    tuple val(meta), path(pytor_files)
    val bin_sizes
    val output_format

    output:
    tuple val(meta), path("*.vcf"), emit: vcf      , optional: true
    tuple val(meta), path("*.tsv"), emit: tsv      , optional: true
    tuple val(meta), path("*.xls"), emit: xls      , optional: true
    tuple val("${task.process}"), val('cnvpytor'), val('1.3.2'), emit: versions_cnvpytor, topic: versions
    // cnvpytor version is hardcoded due to this error when calling cnvpytor --version
    // > cnvpytor --version
    // 2026-03-13 16:14:08,088 - cnvpytor - ERROR - Some reference genome resource files are missing.
    // Run 'cnvpytor -download' as same user who has installed cnvpytor.

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
