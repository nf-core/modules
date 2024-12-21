process AGAT_SPADDINTRONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.0--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config_param = config ? "--config $config" : ""
    def prefix = meta.id ?: gff.getBaseName()
    output = "${prefix}.intron.gff"
    """
    agat_sp_add_introns.pl \\
        --gff $gff \\
        $config_param \\
        --out $output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_add_introns.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = meta.id ?: gff.getBaseName()
    output = "${prefix}.intron.gff"
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_add_introns.pl --help | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
