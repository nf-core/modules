process GT_GFF3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'quay.io/biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gt.gff3")  , emit: gt_gff3  , optional: true
    tuple val(meta), path("*.error.log"), emit: error_log, optional: true
    tuple val("${task.process}"), val('genometools'), eval("gt --version | sed '1!d;s/.* //'"), emit: versions_gt, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gt \\
        gff3 \\
        ${args} \\
        ${gff3} \\
        > "${prefix}.gt.gff3" \\
        2>| >(tee "${prefix}.error.log" >&2) \\
        || echo "Errors from gt-gff3 printed to ${prefix}.error.log"

    if grep -q "gt gff3: error:" "${prefix}.error.log"; then
        echo "gt-gff3 failed to parse ${gff3}"

        rm ${prefix}.gt.gff3
    else
        echo "gt-gff3 successfully parsed ${gff3}"

        mv \\
            ${prefix}.error.log \\
            gt_gff3.stderr
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gt.gff3
    touch ${prefix}.error.log
    touch gt_gff3.stderr
    """
}
