process PASTY {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pasty=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pasty:1.0.0--hdfd78af_0':
        'biocontainers/pasty:1.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv")        , emit: tsv
    tuple val(meta), path("${prefix}.blastn.tsv") , emit: blast
    tuple val(meta), path("${prefix}.details.tsv"), emit: details
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pasty \\
        $args \\
        --prefix $prefix \\
        --assembly $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasty: \$(echo \$(pasty --version 2>&1) | sed 's/^.*pasty, version //;' )
    END_VERSIONS
    """
}
