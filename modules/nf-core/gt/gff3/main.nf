process GT_GFF3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gt.gff3")  , emit: gt_gff3     , optional: true
    tuple val(meta), path("*.error.log"), emit: error_log   , optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gt \\
        gff3 \\
        $args \\
        "$gff3" \\
        > "${prefix}.gt.gff3" \\
        2> >(tee "${prefix}.error.log" >&2) \\
        || echo "Errors from gt-gff3 printed to ${prefix}.error.log"

    if grep -q "gt gff3: error:" "${prefix}.error.log"; then
        echo "gt-gff3 failed to parse $gff3"

        rm \\
            "${prefix}.gt.gff3"
    else
        echo "gt-gff3 successfully parsed $gff3"

        mv \\
            "${prefix}.error.log" \\
            gt_gff3.stderr
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """
}
