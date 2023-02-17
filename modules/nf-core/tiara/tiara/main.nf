process TIARA_TIARA {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "conda-forge::tiara=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiara:1.0.3' :
        'quay.io/biocontainers/tiara:1.0.3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("classification_*.out.txt") , emit: classifications
    tuple val(meta), path("log_*.out.txt")            , emit: log
    tuple val(meta), path("*.fasta*")                 , emit: fasta, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tiara -i ${fasta} \
        -o ${prefix}.out.txt \
        --threads ${task.cpus} \
        ${args}

    mv ${prefix}.out.txt classification_${prefix}.out.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """
}
