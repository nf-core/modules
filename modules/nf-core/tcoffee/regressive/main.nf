process TCOFFEE_REGRESSIVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/t-coffee_pigz:91ac7e26b23bb246':
        'community.wave.seqera.io/library/t-coffee_pigz:7d1373a24f76afe6' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    tuple val(meta3), path(template), path(accessory_information)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def tree_args     = tree ? "-reg_tree $tree" : ""
    def template_args = template ? "-template_file $template" : ""
    def outfile       = compress ? "stdout" : "${prefix}.aln"
    """
    export TEMP='./'
    export TMP_4_TCOFFEE="./"
    export HOME="./"

    t_coffee -reg \
        -seq ${fasta} \
        $tree_args \
        $template_args \
        $args \
        -reg_thread ${task.cpus} \
        -outfile $outfile

    if [ "$compress" = true ]; then
        pigz -cp ${task.cpus} < stdout > ${prefix}.aln.gz
        rm stdout
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    export TMP_4_TCOFFEE="./"
    export HOME="./"

    touch ${prefix}.aln${compress ? '.gz':''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
