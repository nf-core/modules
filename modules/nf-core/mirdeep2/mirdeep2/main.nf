process MIRDEEP2_MIRDEEP2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.2--0':
        'quay.io/biocontainers/mirdeep2:2.0.1.2--0' }"

    input:
    tuple val(meta), path(processed_reads), path(genome_mappings)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(mature), path(hairpin), path(mature_other_species)

    output:
    tuple val(meta), path("result*.{bed,csv,html}")    , emit: outputs
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('mirdeep2'), val("2.0.1"), emit: versions_mirdeep2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mature_species  = mature                ? "${mature}"              : "none"
    def mature_other    = mature_other_species  ? "${mature_other_species}": "none"
    def precursors      = hairpin               ? "${hairpin}"             : "none"

    """
    miRDeep2.pl \\
        ${processed_reads} \\
        ${fasta} \\
        ${genome_mappings} \\
        ${mature_species} \\
        ${mature_other} \\
        ${precursors} \\
        ${args}

    mv result_*.bed result_${prefix}.bed
    mv result_*.csv result_${prefix}.csv
    mv result_*.html result_${prefix}.html
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch result_${prefix}.html
    touch result_${prefix}.bed
    touch result_${prefix}.csv
    """
}
