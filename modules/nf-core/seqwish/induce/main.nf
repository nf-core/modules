process SEQWISH_INDUCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqwish:0.7.11--h5ca1c30_1' :
        'quay.io/biocontainers/seqwish:0.7.11--h5ca1c30_1' }"

    input:
    tuple val(meta), path(paf), path(fasta)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    tuple val("${task.process}"), val("seqwish"), eval("seqwish --version |& sed 's/v//;s/-.*//'"), topic: versions, emit: versions_seqwish

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = paf.join(',') // this ensures that we can actually input a
        // comma-separated list of PAF files as required by
        // https://github.com/nf-core/pangenome. If one wants to use this,
        // ensure that you put a ".collect()" behind your channel.
        // See https://github.com/nf-core/pangenome/blob/34149c6cdc19bce3a7b99f97c769d8986a8d429b/main.nf#L543
        // for an example.
    """
    seqwish \\
        --threads ${task.cpus} \\
        --paf-alns=${input} \\
        --seqs=${fasta} \\
        --gfa=${prefix}.gfa \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gfa
    """

}
