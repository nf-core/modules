process AGAT_SPFILTERFEATUREFROMKILLLIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.0--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    path kill_list
    path config

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def config_param    = config ? "--config $config" : ''
    if( "$gff" == "${prefix}.gff" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    agat_sp_filter_feature_from_kill_list.pl \\
        --gff $gff \\
        --kill_list $kill_list \\
        $config_param \\
        $args \\
        --output "${prefix}.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_filter_feature_from_kill_list.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    if( "$gff" == "${prefix}.gff" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch "${prefix}.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_filter_feature_from_kill_list.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
