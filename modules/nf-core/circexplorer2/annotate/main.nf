process CIRCEXPLORER2_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::circexplorer2=2.3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circexplorer2:2.3.8--pyh864c0ab_1':
        'biocontainers/circexplorer2:2.3.8--pyh864c0ab_1' }"

    input:
    tuple val(meta), path(junctions)
    path(fasta)
    path(gene_annotation)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRCexplorer2 \\
        annotate \\
        -r $gene_annotation \\
        -g $fasta \\
        -b $junctions \\
        -o ${prefix}.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circexplorer2: \$(echo \$(CIRCexplorer2 --version 2>&1) )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circexplorer2: \$(echo \$(CIRCexplorer2 --version 2>&1) )
    END_VERSIONS
    """
}
