process SNAKEMAKE {
    tag "$meta.id"
    label 'process_low'

    // You will have to add all modules to this Conda definition and
    // replace the container definition for one that suits your needs
    conda "bioconda::snakemake=7.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snakemake:7.31.0--hdfd78af_1' :
        'biocontainers/snakemake:7.31.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(inputs)
    tuple val(meta2), path(snakefile)

    output:
    tuple val(meta), path("[!.snakemake|versions.yml]**")         , emit: outputs      , optional: true
    tuple val(meta), path(".snakemake", type: 'dir', hidden: true), emit: snakemake_dir
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cores  = task.cpus ? "--cores ${task.cpus}" : "--cores all"
    """
    snakemake \\
        $args \\
        $cores \\
        --snakefile $snakefile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snakemake: \$(snakemake --version)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cores  = task.cpus ? "--cores ${task.cpus}" : "--cores all"
    """
    snakemake \\
        $args \\
        --snakefile $snakefile \\
        $cores \\
        --dry-run

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snakemake: \$(snakemake --version)
    END_VERSIONS
    """
}
