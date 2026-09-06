process PROOVFRAME_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/proovframe:0.9.7--hdfd78af_1':
        'quay.io/biocontainers/proovframe:0.9.7--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(faa)
    tuple val(meta3), path(db)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('proovframe'), eval("proovframe 2>&1 | sed '/proovframe-v/!d;s/proovframe-v//;s/ .*//'"), topic: versions, emit: versions_proovframe

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_type = db ? "--db ${db}" : (faa ? "--aa ${faa}" : "")

    """
    proovframe \\
        map \\
        ${args} \\
        --threads ${task.cpus} \\
        ${db_type} \\
        -o ${prefix}.tsv \\
        ${fasta}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
