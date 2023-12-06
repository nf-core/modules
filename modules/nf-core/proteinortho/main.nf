process PROTEINORTHO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/proteinortho:6.3.0--h70414c8_0':
        'biocontainers/proteinortho:6.3.0--h70414c8_0' }"

    input:
    tuple val(meta), path(fasta_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${meta.id}.proteinortho.tsv")                     , emit: orthologgroups
    tuple val(meta), path("${meta.id}.proteinortho-graph")                   , emit: orthologgraph
    tuple val(meta), path("${meta.id}.blast-graph")                          , emit: blastgraph
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proteinortho \\
        $args \\
        -cpus=$task.cpus \\
        -project=$prefix \\
        $fasta_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proteinortho : \$(echo \$(proteinortho --version 2>&1) )
        \$(echo \$( "diamond:"(diamond version 2>/dev/null | head -n1 | sed 's/^.*version //;' ) || echo "") )
        \$(echo \$( "blast:"(blastp -version 2>/dev/null |head -n1 | sed 's/^.*: //;'  ) || echo "") )
        \$(echo \$( "mmseqs:"(mmseqp version 2>/dev/null ) || echo "") )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.proteinortho.tsv
    touch ${prefix}.proteinortho-graph
    touch ${prefix}.blast-graph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proteinortho : \$(echo \$(proteinortho --version 2>&1) )
        \$(echo \$( "diamond:"(diamond version 2>/dev/null | head -n1 | sed 's/^.*version //;' ) || echo "") )
        \$(echo \$( "blast:"(blastp -version 2>/dev/null |head -n1 | sed 's/^.*: //;'  ) || echo "") )
        \$(echo \$( "mmseqs:"(mmseqp version 2>/dev/null ) || echo "") )
    END_VERSIONS
    """
}
