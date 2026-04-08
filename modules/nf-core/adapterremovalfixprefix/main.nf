process ADAPTERREMOVALFIXPREFIX {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/adapterremovalfixprefix:0.0.5--hdfd78af_2':
        'biocontainers/adapterremovalfixprefix:0.0.5--hdfd78af_2' }"

    input:
    tuple val(meta), path(fastq)
    output:
    tuple val(meta), path("*.fq.gz"), emit: fixed_fastq
    tuple val("${task.process}"), val('adapterremovalfixprefix'), eval('echo 0.0.5'), emit: versions_adapterremovalfixprefix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${fastq}" == "${prefix}.fq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    AdapterRemovalFixPrefix \\
        ${fastq} \\
        ${args} \\
        | gzip > ${prefix}.fq.gz
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo | gzip > ${prefix}.fq.gz
    """
}
