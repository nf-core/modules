process SNAKEMAKE {
    tag "${meta.id}"
    label 'process_low'

    // You will have to add all modules to this Conda definition and
    // replace the container definition for one that suits your needs
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b91b2eddb4c4c0a5e17721d29764e59a50035b9fa9996eb9cb392829f2d7a533/data'
        : 'community.wave.seqera.io/library/snakemake:9.14.0--dfee75b6201d25c6'}"

    input:
    tuple val(meta), path(inputs)
    tuple val(meta2), path(snakefile)

    output:
    tuple val(meta), path("[!.snakemake|versions.yml]**"), emit: outputs, optional: true
    tuple val(meta), path(".snakemake", type: 'dir', hidden: true), emit: snakemake_dir
    tuple val("${task.process}"), val("snakemake"), eval('snakemake --version'), topic: versions, emit: versions_snakemake

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cores = task.cpus ? "--cores ${task.cpus}" : "--cores all"
    """
    export XDG_CACHE_HOME=\$PWD/snakemake_cache
    mkdir -p \$XDG_CACHE_HOME

    snakemake \\
        ${args} \\
        ${cores} \\
        --snakefile ${snakefile}
    """


    stub:
    def args = task.ext.args ?: ''
    def cores = task.cpus ? "--cores ${task.cpus}" : "--cores all"
    """
    export XDG_CACHE_HOME=\$PWD/snakemake_cache
    mkdir -p \$XDG_CACHE_HOME

    snakemake \\
        ${args} \\
        --snakefile ${snakefile} \\
        ${cores} \\
        --dry-run
    """
}
