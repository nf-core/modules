process GCTA_GSMR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0':
        'quay.io/biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:

    tuple val(meta), path(exposure)
    tuple val(meta2), path(outcome)
    path(reference)

    output:
    path "${meta.id}_${meta2.id}.log", emit: log
    path "${meta.id}_${meta2.id}.gsmr", emit: gsmr
    path "${meta.id}_${meta2.id}.eff_plot.gz", emit: eff_plot, optional: true
    path "${meta.id}_${meta2.id}.mono.badsnps", emit: mono_badsnps, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "${meta.id} ${exposure}" > ${meta.id}.input.txt
    echo "${meta2.id} ${outcome}" > outcome.txt
    file=\$(ls $reference | sed 's/\\.[^.]*\$//')
    echo "${reference}/\$file" | head -n1 > reference.txt

    gcta  \
    --mbfile reference.txt  \
    --gsmr-file ${meta.id}.input.txt outcome.txt \
    --out "${meta.id}_${meta2.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            gcta: \$(gcta 2>&1 | awk '/no analysis has been launched/ {exit 0} {print}' | sed -n 's/.*version \\(v[0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${meta.id}_${meta2.id}.log
    touch ${meta.id}_${meta2.id}.gsmr
    touch ${meta.id}_${meta2.id}.mono.badsnps
    echo "" | gzip > ${meta.id}_${meta2.id}.eff_plot.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            gcta: \$(gcta 2>&1 | awk '/no analysis has been launched/ {exit 0} {print}' | sed -n 's/.*version \\(v[0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
