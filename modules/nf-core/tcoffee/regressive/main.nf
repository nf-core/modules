// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TCOFFEE_REGRESSIVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e85c972da17c3e7d5e8b283cf61b106111180e3fc02d4b04f0933106b4ca50b3/data'
        : 'community.wave.seqera.io/library/t-coffee_pigz:aa6b4f1ff4adb194'}"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    tuple val(meta3), path(template), path(accessory_information)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def tree_args          = tree ? "-reg_tree $tree" : ""
    def template_args      = template ? "-template_file $template" : ""
    def default_out_format = ("-output" in "${args}") ? "" : "-output fasta_aln"
    def outfile            = compress ? "stdout" : "${prefix}.aln"
    """
    export TEMP='./'
    t_coffee -reg \
        -seq ${fasta} \
        $tree_args \
        $template_args \
        $args \
        $default_out_format \
        -reg_thread ${task.cpus} \
        -outfile $outfile

    if [ "$compress" = true ]; then
        pigz -cp ${task.cpus} < stdout > ${prefix}.aln.gz
        rm stdout
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch ${prefix}.aln${compress ? '.gz':''}
    """
}
