process GT_GFF3VALIDATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'quay.io/biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path('*.success.log')  , emit: success_log , optional: true
    tuple val(meta), path('*.error.log')    , emit: error_log   , optional: true
    tuple val("${task.process}"), val('genometools'), eval("gt --version | sed '1!d;s/.* //'"), emit: versions_gt, topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gt \\
        gff3validator \\
        ${gff3} \\
        > ${prefix}.stdout \\
        2>| >(tee "${prefix}.stderr" >&2) \\
        || echo "Errors from gt-gff3validator printed to ${prefix}.error.log"

    if grep -q "input is valid GFF3" "${prefix}.stdout"; then
        echo "Validation successful..."
        # emit stdout to the success output channel
        mv \\
            ${prefix}.stdout \\
            ${prefix}.success.log
    else
        echo "Validation failed..."
        # emit stderr to the error output channel
        mv \\
            ${prefix}.stderr \\
            ${prefix}.error.log
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.success.log
    touch ${prefix}.error.log
    """
}
