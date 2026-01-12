process S4PRED_RUNMODEL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/s4pred:1.2.1--pyhdfd78af_1':
        'biocontainers/s4pred:1.2.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}"), emit: preds
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir ${prefix}

    run_model \\
        $args \\
        --threads $task.cpus \\
        --save-files \\
        --outdir ${prefix} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        s4pred: $VERSION
    END_VERSIONS
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p ${prefix}
    touch ${prefix}/test.ss2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        s4pred: $VERSION
    END_VERSIONS
    """
}
