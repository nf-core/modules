process SCOARY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::scoary=1.6.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scoary:1.6.16--py_2' :
        'quay.io/biocontainers/scoary:1.6.16--py_2' }"

    input:
    tuple val(meta), path(genes), path(traits)
    path(tree)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def newick_tree = tree ? "-n ${tree}" : ""
    """
    scoary \\
        $args \\
        --no-time \\
        --threads $task.cpus \\
        --traits $traits \\
        --genes $genes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scoary: \$( scoary --version 2>&1 )
    END_VERSIONS
    """
}
