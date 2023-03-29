
process LOFREQ3_PILEUP {
    tag "$meta.id"
    label 'process_low'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/vojalu/lofreq:3.0' }"

    input:
    tuple val(meta),
          path(lofreq3_dir)

    output:
    tuple val(meta), path("*.json"), emit: pileup
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    lofreq call -p

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq3: \$(echo "3.0 - reimplementation of 2 in Nim")
    END_VERSIONS
    """
}
