process MODKIT_EXTRACT_CALLS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tsv{,.gz}"), emit: tsv
    tuple val(meta), path("*.log")      , emit: log, optional: true
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def reference  = fasta ? "--reference ${fasta}" : ''
    def out_suffix = args.tokenize().contains('--bgzf') ? 'tsv.gz' : 'tsv'
    """
    modkit \\
        extract \\
        calls \\
        $args \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${bam} \\
        ${prefix}.${out_suffix}
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def out_suffix = args.tokenize().contains('--bgzf') ? 'tsv.gz' : 'tsv'
    """
    touch ${prefix}.${out_suffix}
    """
}
