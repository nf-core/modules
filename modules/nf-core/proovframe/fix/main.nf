process PROOVFRAME_FIX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/proovframe:0.9.7--hdfd78af_1':
        'biocontainers/proovframe:0.9.7--hdfd78af_1' }"

    input:
    tuple val(meta) , path(fa)
    tuple val(meta2), path(tsv)

    output:
    tuple val(meta), path("*.fa"), emit: out_fa
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proovframe   \\
        fix \\
        ${args} \\
        -o ${prefix}.fa  \\
        ${fa} \\
        ${tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proovframe: \$(proovframe 2>&1 | grep -o 'proovframe-v[0-9]*\\.[0-9]*\\.[0-9]*' | grep -o '[0-9]*\\.[0-9]*\\.[0-9]*')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proovframe: \$(proovframe 2>&1 | grep -o 'proovframe-v[0-9]*\\.[0-9]*\\.[0-9]*' | grep -o '[0-9]*\\.[0-9]*\\.[0-9]*')
    END_VERSIONS
    """
}
