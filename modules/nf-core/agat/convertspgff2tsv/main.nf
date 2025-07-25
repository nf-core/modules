process AGAT_CONVERTSPGFF2TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_convert_sp_gff2tsv.pl \\
        --gff ${gff} \\
        --output ${prefix}.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gff2tsv.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gff2tsv.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
