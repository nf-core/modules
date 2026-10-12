process GENIN2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genin2:2.1.6--pyhdfd78af_0':
        'quay.io/biocontainers/genin2:2.1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('genin2'), eval("genin2 --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_genin2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genin2 \\
        $args \\
        -o ${prefix}.tsv \\
        $fasta
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tsv
    """
}
