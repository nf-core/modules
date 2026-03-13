process TCOFFEE_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/75830929ed209853ab621dafcc4105618eba02408a731e4fd4aeee94ada424bc/data':
        'community.wave.seqera.io/library/t-coffee_tmalign_pigz:be7dac2ae6aba380' }"

    input:
    tuple val(meta) ,  path(fasta)
    tuple val(meta2),  path(tree)
    tuple val(meta3),  path(template), path(accessory_information)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    // in the args there might be the request to generate a lib file, so the following is an optional output
    tuple val(meta), path("*.*lib")     , emit: lib, optional : true

    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | awk \'{gsub("Version_", ""); print \\$3}\''), emit: versions_tcoffee, topic: versions
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
