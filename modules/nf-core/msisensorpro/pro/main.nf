process MSISENSORPRO_PRO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro%3A1.3.0--hfef96ef_0':
        'quay.io/biocontainers/msisensor-pro:1.3.0--hfef96ef_0' }"

    input:
    tuple val(meta), path(cram)
    path("${cram}.crai")
    path(list)
    path(fasta)
    path("${fasta}.fai")

    output:
    tuple val(meta), path("${prefix}")          , emit: summary_msi
    tuple val(meta), path("${prefix}_all")      , emit: all_msi
    tuple val(meta), path("${prefix}_dis")      , emit: dis_msi
    tuple val(meta), path("${prefix}_unstable") , emit: unstable_msi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        pro \\
        -d $list \\
        -t $cram \\
        -g $fasta \\
        -o ${prefix} \\
        -b $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensorpro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\s*//p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix} ${prefix}_all ${prefix}_dis ${prefix}_unstable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensorpro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\s*//p')
    END_VERSIONS
    """
}
