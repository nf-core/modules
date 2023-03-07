process COOLER_CLOAD {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cooler=0.8.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(pairs), path(index), val(cool_bin)
    path chromsizes

    output:
    tuple val(meta), path("*.cool"), val(cool_bin), emit: cool
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nproc  = args.contains('pairix') || args.contains('tabix')? "--nproc $task.cpus" : ''

    """
    cooler cload \\
        $args \\
        $nproc \\
        ${chromsizes}:${cool_bin} \\
        $pairs \\
        ${prefix}.cool

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
