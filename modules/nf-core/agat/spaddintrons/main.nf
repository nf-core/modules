process AGAT_SPADDINTRONS {
    tag '$gff'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.2.0--pl5321hdfd78af_0':
        'biocontainers/agat:1.2.0--pl5321hdfd78af_0' }"

    input:
    path gff
    path config

    output:
    path "${output}"   , emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config_param = config ? "--config $config" : ""
    output = "${gff.getBaseName()}.intron.gff"
    """
    agat_sp_add_introns.pl \\
        --gff $gff \\
        $config_param \\
        --out $output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_add_introns.pl --help | sed '3!d; s/.*v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    output = "output.intron.gff"
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_add_introns.pl --help | sed '3!d; s/.*v//')
    END_VERSIONS
    """
}
