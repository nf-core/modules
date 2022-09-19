process SOMALIER_RELATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::somalier=0.2.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/somalier:0.2.15--h37c5b7d_0':
        'quay.io/biocontainers/somalier:0.2.15--h37c5b7d_0' }"

    input:
    tuple val(meta), path(extract)
    path(sample_groups)
    path(ped)

    output:
    tuple val(meta), path("*.html"),          emit: html
    tuple val(meta), path("*.pairs.tsv"),     emit: pairs_tsv
    tuple val(meta), path("*.samples.tsv"),   emit: samples_tsv
<<<<<<< HEAD
    path "versions.yml",                                  emit: versions
=======
    path "versions.yml",                      emit: versions
>>>>>>> 6c915826 (requested changes)

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
<<<<<<< HEAD
    -o ${prefix} \\
    ${input_list} \\
    ${args} \\
    ${sample_groups_command} \\
    ${ped_command}
=======
        -o ${prefix} \\
        ${input_list} \\
        ${args} \\
        ${sample_groups_command} \\
        ${ped_command}
>>>>>>> 6c915826 (requested changes)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """

}
