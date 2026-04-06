process ANNOSINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annosine2:2.0.8--pyh7e72e81_0':
        'biocontainers/annosine2:2.0.8--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(fasta)
    val mode

    output:
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val(meta), path("${prefix}.fa") , emit: fa, optional: true
    tuple val("${task.process}"), val('annosine'), eval("pip show annosine2 | sed -n 's/Version: //p'"), emit: versions_annosine, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}_annosine"
    if ( "${fasta}" == "${prefix}.fa" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    AnnoSINE_v2 \\
        ${args} \\
        --threads ${task.cpus} \\
        ${mode} \\
        ${fasta} \\
        ${prefix} \\
        &> >(tee ${prefix}.log 2>&1)

    mv \\
        ${prefix}/Seed_SINE.fa \\
        ${prefix}.fa \\
        || echo 'AnnoSINE_v2 did not find SINE sequences. See log for details!'
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}_annosine"
    if ( "${fasta}" == "${prefix}.fa" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.log
    touch ${prefix}.fa
    """
}
