// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TCOFFEE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed121be7a85e6c4f964ffcda6a15ce7f341bbe6050b921e39cb69122ed34ed7b/data'
        : 'community.wave.seqera.io/library/t-coffee_tmalign_pigz:c536fc45adb585e9'}"

    input:
    tuple val(meta) ,  path(fasta)
    tuple val(meta2),  path(tree)
    tuple val(meta3),  path(template), path(accessory_information)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    // in the args there might be the request to generate a lib file, so the following is an optional output
    tuple val(meta), path("*.*lib")     , emit: lib, optional : true

    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def tree_args          = tree ? "-usetree $tree" : ""
    def template_args      = template ? "-template_file $template" : ""
    def outfile            = compress ? "stdout" : "${prefix}.aln"
    def default_out_format = ("-output" in "${args}") ? "" : "-output fasta_aln"
    def write_output       = compress ? " | pigz -cp ${task.cpus} > ${prefix}.aln.gz" : ""
    """
    export TEMP='./'
    export TMP_4_TCOFFEE="./"
    export HOME="./"
    t_coffee -seq ${fasta} \
        $tree_args \
        $template_args \
        $args \
        $default_out_format \
        -thread ${task.cpus} \
        -outfile $outfile \
        $write_output
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    export TMP_4_TCOFFEE="./"
    export HOME="./"
    touch ${prefix}.aln${compress ? '.gz':''}
    """
}
