process MIRDEEP2_MIRDEEP2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.2--0':
        'biocontainers/mirdeep2:2.0.1.2--0' }"

    input:
    tuple val(meta), path(processed_reads), path(genome_mappings)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(mature), path(hairpin), path(mature_other_species)

    output:
    tuple val(meta), path("result*.{bed,csv,html}")    , emit: outputs
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.0.1'
    def mature_species  = mature                ? "${mature}"              : "none"
    def mature_other    = mature_other_species  ? "${mature_other_species}": "none"
    def precursors      = hairpin               ? "${hairpin}"             : "none"

    """
    miRDeep2.pl \\
        $processed_reads \\
        $fasta \\
        $genome_mappings \\
        $mature_species \\
        $mature_other \\
        $precursors \\
        $args

    mv result_*.bed result_${prefix}.bed
    mv result_*.csv result_${prefix}.csv
    mv result_*.html result_${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirdeep2: \$(echo "$VERSION")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.0.1'
    """
    touch result_${prefix}.html
    touch result_${prefix}.bed
    touch result_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirdeep2: \$(echo "$VERSION")
    END_VERSIONS
    """
}
