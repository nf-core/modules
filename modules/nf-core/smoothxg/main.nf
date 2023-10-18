process SMOOTHXG {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::smoothxg=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoothxg:0.7.1--h40c17d1_0' :
        'biocontainers/smoothxg:0.7.1--h40c17d1_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*smoothxg.gfa"), emit: gfa
    path("*.maf") , optional: true, emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    smoothxg \\
        --threads=$task.cpus \\
        --gfa-in=${gfa} \\
        --smoothed-out=${prefix}.smoothxg.gfa \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smoothxg: \$(smoothxg --version 2>&1 | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
