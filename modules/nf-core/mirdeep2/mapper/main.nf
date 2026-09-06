process MIRDEEP2_MAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.2--0':
        'quay.io/biocontainers/mirdeep2:2.0.1.2--0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index, stageAs: '*')

    output:
    tuple val(meta), path('*.fa'), path('*.arf'), emit: outputs
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('mirdeep2'), val("2.0.1"), emit: versions_mirdeep2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mapper.pl \\
        ${reads} \\
        $args \\
        -p ${index}/${meta2.id}  \\
        -s ${prefix}_collapsed.fa \\
        -t ${prefix}_reads_collapsed_vs_${meta2.id}_genome.arf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}reads_vs_refdb.arf
    """
}
