process TCOFFEE_REGRESSIVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
        'biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    tuple val(meta3), path(template), path(accessory_information)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | awk \'{gsub("Version_", ""); print \\$3}\''), emit: versions_tcoffee, topic: versions
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
