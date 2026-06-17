process VCLUST_ALIGN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'quay.io/biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(filter)
    val save_alignment

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("*.ids.tsv")    , emit: ids
    tuple val(meta), path("*.aln.tsv")    , emit: aln, optional: true
    tuple val("${task.process}"), val('vclust'), eval("vclust --version | sed 's/v//'"), topic: versions, emit: versions_vclust


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def filter_argument = filter ? "--filter ${filter}" : ''
    save_alignment_option = save_alignment ? "--out-aln ${prefix}.aln.tsv" : ''
    """
    vclust \\
        align \\
        $args \\
        ${filter_argument} \\
        ${save_alignment_option} \\
        -t $task.cpus \\
        -i ${fasta} \\
        -o ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.ids.tsv
    touch ${prefix}.aln.tsv
    """
}
