process JVARKIT_WGSCOVERAGEPLOTTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    output:
    tuple val(meta),  path("*.svg"), emit: output
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p TMP

    jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP wgscoverageplotter \\
        -R ${fasta} \\
        ${args} \\
        ${bam} > "${prefix}.svg"

    rm -rf TMP

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.svg"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """
}
