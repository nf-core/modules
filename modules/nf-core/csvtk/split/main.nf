process CSVTK_SPLIT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/917edb71b915f07fa2838c20e3c731181d3d315cbf8a9bfead41412d2b4ae062/data' :
        'community.wave.seqera.io/library/csvtk:0.37.0--113625988dd3285d' }"

    input:
    tuple val(meta), path(csv)
    val in_format
    val out_format

    output:
    tuple val(meta), path("*.${out_extension}"), emit: split_csv
    tuple val("${task.process}"), val('csvtk'), eval("csvtk version | sed -e 's/csvtk v//g'"), emit: versions_csvtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def delimiter = in_format == "tsv" ? "--tabs" : (in_format == "csv" ? "--delimiter ',' " : in_format)
    def out_delimiter = out_format == "tsv" ? "--out-tabs" : (out_format == "csv" ? "--out-delimiter ',' " : out_format)
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    """
    sed -i.bak '/^##/d' $csv
    csvtk \\
        split \\
        ${args} \\
        --num-cpus ${task.cpus} \\
        --out-prefix ${prefix} \\
        ${delimiter} \\
        ${out_delimiter} \\
        ${csv}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = args.contains('--out-delimiter "\t"') || args.contains('-D "\t"') || args.contains("-D \$'\t'") ? "tsv" : "csv"
    """
    touch ${prefix}.${out_extension}
    """
}
