process PROOVFRAME_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/proovframe:0.9.7--hdfd78af_1':
        'biocontainers/proovframe:0.9.7--hdfd78af_1' }"

    input:
    tuple val(meta), path(faa), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proovframe \\
        map \\
        ${args} \\
        -a ${faa} \\
        -o ${prefix}.tsv  \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proovframe: \$(proovframe 2>&1 | grep -o 'proovframe-v[0-9]*\\.[0-9]*\\.[0-9]*' | grep -o '[0-9]*\\.[0-9]*\\.[0-9]*')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proovframe: \$(proovframe 2>&1 | grep -o 'proovframe-v[0-9]*\\.[0-9]*\\.[0-9]*' | grep -o '[0-9]*\\.[0-9]*\\.[0-9]*')
    END_VERSIONS
    """
}
