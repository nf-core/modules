process MAPAD_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapad:0.44.1--ha96b9cd_0':
        'biocontainers/mapad:0.44.1--ha96b9cd_0' }"

    input:
    tuple val(meta), path(fasta, stageAs: "mapad/*")

    output:
    tuple val(meta), path("mapad/"), emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mapad \\
        index \\
        $args \\
        --reference $fasta \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapad: \$(echo \$(mapad --version) | sed 's/^mapAD //' )
    END_VERSIONS
    """
}
