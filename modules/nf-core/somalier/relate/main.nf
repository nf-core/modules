
process SOMALIER_RELATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/somalier:0.2.19--h0c29559_0':
        'biocontainers/somalier:0.2.19--h0c29559_0' }"

    input:
    tuple val(meta), path(extract), path(ped)
    path(sample_groups)

    output:
    tuple val(meta), path("*.html"),          emit: html
    tuple val(meta), path("*.pairs.tsv"),     emit: pairs_tsv
    tuple val(meta), path("*.samples.tsv"),   emit: samples_tsv
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_list = extract.collect{"$it"}.join(' ')
    def prefix = task.ext.prefix ?: "$meta.id"
    def sample_groups_command = sample_groups ? "-g $sample_groups" : ""
    def ped_command = ped ? "-p $ped" : ""

    """
    somalier relate \\
        -o ${prefix} \\
        ${input_list} \\
        ${args} \\
        ${sample_groups_command} \\
        ${ped_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"

    """
    touch ${prefix}.html
    touch ${prefix}.pairs.tsv
    touch ${prefix}.samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """
}
