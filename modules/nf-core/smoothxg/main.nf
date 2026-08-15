process SMOOTHXG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoothxg:0.8.0--h40c17d1_0' :
        'quay.io/biocontainers/smoothxg:0.8.0--h40c17d1_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*smoothxg.gfa"), emit: gfa
    path("*.maf") , optional: true, emit: maf
    tuple val("${task.process}"), val('smoothxg'), eval("smoothxg --version 2>&1 | sed 's/^v//; s/-.*//'"  ), topic: versions, emit: versions_smoothxg

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.smoothxg.gfa
    """
}
