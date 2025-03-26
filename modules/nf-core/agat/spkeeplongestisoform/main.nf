process AGAT_SPKEEPLONGESTISOFORM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gxf)
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
    output = "${prefix}.longest.gff"
    """
    agat_sp_keep_longest_isoform.pl \\
        --gff $gxf \\
        $config_param \\
        --out $output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat --version)
    END_VERSIONS
    """

    stub:
    def prefix = meta.id ?: gff.getBaseName()
    output = "${prefix}.longest.gff"
    """
    touch ${output}
    touch ${gxf}.agat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat --version)
    END_VERSIONS
    """
}