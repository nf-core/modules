process UPP_ALIGN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sepp:4.5.6--py312h87e0c26_4'
        : 'biocontainers/sepp:4.5.6--py312h87e0c26_4'}"

    input:
    tuple val(meta), path(fasta_unaligned), path(fasta_aligned)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    tuple val("${task.process}"), val('sepp'), eval('run_upp.py -v | grep "run_upp" | cut -f2 -d" "'), emit: versions_sepp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tree_args = tree ? "-t ${tree}" : ""
    def seq_cmd = fasta_unaligned ? "--sequence_file ${fasta_unaligned}" : ""
    def align_cmd = fasta_aligned ? "--alignment ${fasta_aligned}" : ""

    """
    if [ "${workflow.containerEngine}" = 'singularity' ]; then
        export CONDA_PREFIX="/usr/local/"
        export PASTA_TOOLS_DEVDIR="/usr/local/bin/"
    fi

    run_upp.py \\
        ${args} \\
        ${tree_args} \\
        ${seq_cmd} \\
        ${align_cmd} \\
        -x ${task.cpus} \\
        -p ./upp_tmp/ \\
        -d ./upp_output/ \\
        -o ${prefix}

    mv upp_output/${prefix}_alignment.fasta ${prefix}.aln

    if ${compress}; then
        gzip ${prefix}.aln
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ "${compress}" = true ]; then
        echo | gzip > "${prefix}.aln.gz"
    else
        touch "${prefix}.aln"
    fi
    """
}
