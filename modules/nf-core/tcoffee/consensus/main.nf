process TCOFFEE_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/t-coffee_pigz:f47b85d70360f1a0':
        'community.wave.seqera.io/library/t-coffee_pigz:6c9b2f8b97ee55e5' }"


    input:
    tuple val(meta) , path(aln)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.{aln,aln.gz}"), emit: alignment
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tree_args = tree ? "-usetree $tree" : ""
    def outfile = compress ? "stdout" : "${prefix}.aln"
    def write_output = compress ? " | pigz -cp ${task.cpus} > ${prefix}.aln.gz" : ""
    """
    export TEMP='./'
    t_coffee -aln ${aln} \
        $tree_args \
        $args \
        -thread ${task.cpus} \
        -outfile $outfile \
        $write_output

    if [ -f stdout ] && [ "$compress" = true ]; then
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
    export TEMP='./'
    touch ${prefix}.aln${compress ? '.gz':''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
