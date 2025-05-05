process FUSIONCATCHER_FUSIONCATCHER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_5':
        'biocontainers/fusioncatcher:1.33--hdfd78af_5' }"

    input:
    tuple val(meta), path(fastqs)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.fusion-genes.txt"), emit: fusions, optional: true
    tuple val(meta), path("*.summary.txt")     , emit: summary, optional: true
    tuple val(meta), path("*.log")             , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = fastqs.join(",")

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[FusionCatcher] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    fusioncatcher \\
        --input=${input} \\
        --output=. \\
        --data=${reference} \\
        --threads=${task.cpus} \\
        --Xmx=${avail_mem}M \\
        ${args}

    mv final-list_candidate-fusion-genes.txt ${prefix}.fusion-genes.txt
    mv summary_candidate_fusions.txt ${prefix}.summary.txt
    mv fusioncatcher.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fusion-genes.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
