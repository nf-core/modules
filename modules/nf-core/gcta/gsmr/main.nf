process GCTA_GSMR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0':
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta) , path(exposure)
    tuple val(meta2), path(outcome)
    path(reference)

    output:
    tuple val(meta), val(meta2), path("*.log")           , emit: log
    tuple val(meta), val(meta2), path("*.gsmr")          , emit: gsmr
    tuple val(meta), val(meta2), path("*.eff_plot.gz")   , emit: eff_plot, optional: true
    tuple val(meta), val(meta2), path("*.mono.badsnps"), emit: mono_badsnps, optional: true
    tuple val("${task.process}"), val('gcta'), eval('gcta 2>&1 | grep -oE "v[0-9]+\\.[0-9]+\\.[0-9]+" | sed \'s/v//\''), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta2.id}"
    """
    echo "${meta.id} ${exposure}" > ${meta.id}.input.txt
    echo "${meta2.id} ${outcome}" > outcome.txt
    file=\$(ls $reference | sed 's/\\.[^.]*\$//')
    echo "${reference}/\$file" | head -n1 > reference.txt

    gcta  \\
        $args \\
        --mbfile reference.txt  \\
        --gsmr-file ${meta.id}.input.txt outcome.txt \\
        --out "${prefix}"
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta2.id}"
    """
    touch ${prefix}.log
    touch ${prefix}.gsmr
    touch ${prefix}.mono.badsnps
    echo "" | gzip > ${prefix}.eff_plot.gz

    """
}
