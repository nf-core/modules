process DEDUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/dedup:0.12.9--hdfd78af_0'
        : 'quay.io/biocontainers/dedup:0.12.9--hdfd78af_0'}"

    input:
    tuple val(meta), path(bam)

    output:
    // _rmdup is hardcoded output from dedup
    tuple val(meta), path("${prefix}.bam"),  emit: bam
    tuple val(meta), path("${prefix}.json"), emit: json
    tuple val(meta), path("${prefix}.hist"), emit: hist
    tuple val(meta), path("${prefix}.log"),  emit: log
    tuple val("${task.process}"), val('dedup'), eval("dedup --version | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+'"), emit: versions_dedup, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bam_prefix = bam.baseName
    """
    dedup \\
        -Xmx${task.memory.toGiga()}g  \\
        -i ${bam} \\
        -o . \\
        ${args}

    mv ${bam_prefix}_rmdup.bam ${prefix}.bam
    mv *.json ${prefix}.json
    mv *.hist ${prefix}.hist
    mv *.log ${prefix}.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    touch ${prefix}.hist
    touch ${prefix}.log
    touch ${prefix}.bam
    """
}
