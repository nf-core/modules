process UPP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/sepp_pigz:d72591720d0277b1':
        'community.wave.seqera.io/library/sepp_pigz:ea6dbc7704a2e251' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tree_args = tree ? "-t $tree" : ""
    """

    if [ "$workflow.containerEngine" = 'singularity' ]; then
        export CONDA_PREFIX="/opt/conda/"
        export PASTA_TOOLS_DEVDIR="/opt/conda/bin/"
    fi

    run_upp.py \\
        $args \\
        -x $task.cpus \\
        -s ${fasta} \\
        -d . \\
        -o ${prefix} \\
        -p ./upp-temporary

    mv ${prefix}_alignment.fasta ${prefix}.aln

    if ${compress}; then
        pigz -p ${task.cpus} ${prefix}.aln
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upp: \$(run_upp.py -v | grep "run_upp" | cut -f2 -d" ")
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """


    if [ "$compress" = true ]; then
        echo | gzip > "${prefix}.aln.gz"
    else
        touch "${prefix}.aln"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upp: \$(run_upp.py -v | grep "run_upp" | cut -f2 -d" ")
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
