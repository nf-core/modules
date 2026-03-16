process RNAQUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rnaquast:2.3.0--h9ee0642_0':
        'biocontainers/rnaquast:2.3.0--h9ee0642_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    reference = reference ? "--reference ${reference}" : ''
    gtf = gtf ? "--gtf ${gtf}" : ''
    """
    rnaQUAST.py \\
        $args \\
        --threads $task.cpus \\
        --transcripts ${fasta} \\
        ${reference} \\
        ${gtf} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaquast: \$(rnaQUAST.py -h | grep -i 'rnaQUAST.py v' | sed -E 's/.*v\\.([0-9.]+).*/\\1/')
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/rnaQUAST.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaquast: \$(rnaQUAST.py -h | grep -i 'rnaQUAST.py v' | sed -E 's/.*v\\.([0-9.]+).*/\\1/')
    END_VERSIONS
    """
}
