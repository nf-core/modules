process VCONTACT3_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcontact3:3.1.3--pyhdfd78af_0' :
        'biocontainers/vcontact3:3.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(database)
    val is_nucleotide

    output:
    tuple val(meta), path("${prefix}/"), emit: out
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_command = is_nucleotide ? "-n ${fasta}" : "-p ${fasta}"
    """
    vcontact3 \\
        run \\
        $args \\
        -t $task.cpus \\
        ${input_command} \\
        -d ${database} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact3: \$(vcontact3 version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact3: \$(vcontact3 version)
    END_VERSIONS
    """
}
