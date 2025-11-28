process ORFIPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orfipy:0.0.4--py310h184ae93_4':
        'biocontainers/orfipy:0.0.4--py310h184ae93_4' }"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("orfipy_outs/${prefix}.bed"), emit: bed
    path("versions.yml")                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    orfipy \\
        ${infile} \\
        --outdir orfipy_outs \\
        --bed ${prefix}.bed \\
        --procs ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orfipy: \$(orfipy --version | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p orfipy_outs/
    touch orfipy_outs/${prefix}.bed

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orfipy: \$(orfipy --version | sed 's/.*version //')
    END_VERSIONS
    """
}