process GANON_BUILDCUSTOM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ganon=1.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ganon:1.5.1--py310h8abeb55_0':
        'quay.io/biocontainers/ganon:1.5.1--py310h8abeb55_0' }"

    input:
    tuple val(meta), path(input)
    path taxonomy_files
    path genome_size_files

    output:
    tuple val(meta), path("*.{ibf,tax}")          , emit: db
    tuple val(meta), path("*.info.tsv")           , emit: info            , optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxonomy_args     = taxonomy_files    ? "--taxonomy-files ${taxonomy_files}" : ""
    def genome_size_args  = genome_size_files ? "--genome-size-files ${genome_size_files}" : ""
    """
    ganon \\
        build-custom \\
        --threads ${task.cpus} \\
        --input $input \\
        --db-prefix ${prefix} \\
        $taxonomy_args \\
        $genome_size_args \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxonomy_args     = taxonomy_files    ? "--taxonomy-files ${taxonomy_files}" : ""
    def genome_size_args  = genome_size_files ? "--genome-size-files ${genome_size_files}" : ""
    """
    touch ${prefix}.ibf
    touch ${prefix}.tax
    touch ${prefix}.info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """
}
