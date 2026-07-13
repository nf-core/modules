process PRIORCONS_BUILDPRIORS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/priorcons:0.1.0--pyhdfd78af_0' :
    'biocontainers/priorcons:0.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(alignment)  // Input alignment FASTA file
    val   ref_id                       // Reference sequence ID

    output:
    tuple val(meta), path("${meta.id}.parquet"), emit: priors
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    priorcons \\
        build-priors \\
        --input ${alignment} \\
        --ref ${ref_id} \\
        --output ${prefix}.parquet \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        priorcons: \$(priorcons --version 2>&1 | sed 's/priorcons //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        priorcons: \$(priorcons --version 2>&1 | sed 's/priorcons //')
    END_VERSIONS
    """
}
