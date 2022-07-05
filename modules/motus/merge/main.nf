VERSION = '3.0.1'

process MOTUS_MERGE {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::motus=3.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0' }"

    input:
    path input
    path db // to stop docker saying it can't find it... would have to have the module in upstream steps anyway
    path profile_version_yml, stageAs: 'profile_version.yml'
    val biom_format

    output:
    path("*.txt") , optional: true, emit: txt
    path("*.biom"), optional: true, emit: biom
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = 'motus_merged'
    def cmd_input = input.size() > 1 ? "-i ${input.join(',')}" : input.isDirectory() ? "-d ${input}" : "-i ${input}"
    def output = biom_format ? "-B -o ${prefix}.biom" : "-o ${prefix}.txt"
    """
    motus \\
        merge \\
        -db $db \\
        ${cmd_input} \\
        $args \\
        ${output}

    ## Take version from the mOTUs/profile module output, as cannot reconstruct
    ## version without having database staged in this directory.
    VERSION=\$(cat ${profile_version_yml} | grep '/*motus:.*' | sed 's/.*otus: //g')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$VERSION
    END_VERSIONS
    """
}
