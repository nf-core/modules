process MOTUS_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0':
        'biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path db // to stop docker saying it can't find it... would have to have the module in upstream steps anyway
    path profile_version_yml, stageAs: 'profile_version.yml'

    output:
    tuple val(meta), path("*.txt") , optional: true, emit: txt
    tuple val(meta), path("*.biom"), optional: true, emit: biom
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd_input = input.size() > 1 ? "-i ${input.join(',')}" : input.isDirectory() ? "-d ${input}" : "-i ${input}"
    def suffix = task.ext.args?.contains("-B") ? "biom" : "txt"
    """
    motus \\
        merge \\
        -db $db \\
        ${cmd_input} \\
        $args \\
        -o ${prefix}.${suffix}

    ## Take version from the mOTUs/profile module output, as cannot reconstruct
    ## version without having database staged in this directory.
    VERSION=\$(cat ${profile_version_yml} | grep '/*motus:.*' | sed 's/.*otus: //g')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd_input = input.size() > 1 ? "-i ${input.join(',')}" : input.isDirectory() ? "-d ${input}" : "-i ${input}"
    def suffix = task.ext.args?.contains("-B") ? "biom" : "txt"

    """
    touch ${prefix}.txt

    VERSION=\$(cat ${profile_version_yml} | grep '/*motus:.*' | sed 's/.*otus: //g')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$VERSION
    END_VERSIONS
    """

}
