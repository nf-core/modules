process PBPTYPER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pbptyper=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbptyper:1.0.2--hdfd78af_0':
        'biocontainers/pbptyper:1.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("*.tblastn.tsv"), emit: blast
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_args = db ? '--db ${db}' : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbptyper \\
        $db_args \\
        $args \\
        --prefix $prefix \\
        --assembly $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbptyper: \$(echo \$(pbptyper --version 2>&1) | sed 's/^.*pbptyper, version //;' )
    END_VERSIONS
    """
}
