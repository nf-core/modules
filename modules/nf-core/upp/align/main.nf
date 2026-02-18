process UPP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    //TODO: Use bioconda container once the recipe is fixed (see https://github.com/smirarab/sepp/issues/141)
    container "nf-core/multiplesequencealign_upp2:4.5.5"

    input:
    tuple val(meta) , path(fasta_unaligned), path(fasta_aligned)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    tuple val("${task.process}"), val('upp'), eval('run_upp.py -v | grep "run_upp" | cut -f2 -d" "'), emit: versions_upp, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tree_args = tree ? "-t $tree" : ""
    def seq_cmd = fasta_unaligned ? "--sequence_file ${fasta_unaligned}" : ""
    def align_cmd = fasta_aligned ? "--alignment ${fasta_aligned}" : ""

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Upp align module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    """
    if [ "$workflow.containerEngine" = 'singularity' ]; then
        export CONDA_PREFIX="/opt/conda/"
        export PASTA_TOOLS_DEVDIR="/opt/conda/bin/"
    fi

    run_upp.py \\
        ${args} \\
        ${tree_args} \\
        -x ${task.cpus} \\
        ${seq_cmd} \\
        ${align_cmd} \\
        -d . \\
        -o ${prefix} \\
        -p ./upp-temporary

    mv ${prefix}_alignment.fasta ${prefix}.aln

    if ${compress}; then
        pigz -p ${task.cpus} ${prefix}.aln
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
