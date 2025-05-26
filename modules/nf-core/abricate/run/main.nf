process ABRICATE_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/abricate%3A1.0.1--ha8f3691_1'
        : 'biocontainers/abricate:1.0.1--ha8f3691_1'}"

    input:
    tuple val(meta), path(assembly)
    path databasedir

    output:
    tuple val(meta), path("*.txt"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def datadir = databasedir ? "--datadir ${databasedir}" : ''
    if ("${assembly}" == "${prefix}.fasta") {
        error("File name clash in internal file preparation, please modify prefix in module configuration to disambiguate!")
    }
    """
    ## Symlink to rename the file to allow specifying the prefix variable  inside report
    ## As the variable is what is used as the sample ID in the report file
    ln -s ${assembly} ${prefix}.fasta

    abricate \\
        ${prefix}.fasta \\
        ${args} \\
        ${datadir} \\
        --threads ${task.cpus} \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def datadir = databasedir ? '--datadir ${databasedir}' : ''
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
