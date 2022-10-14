process CIRCEXPLORER2_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::circexplorer2=2.3.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circexplorer2:2.3.8--pyh864c0ab_1':
        'quay.io/biocontainers/circexplorer2:2.3.8--pyh864c0ab_1' }"

    input:
    tuple val(meta), path(fusions)

    output:
    tuple val(meta), path("*.bed"), emit: junction
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligner = "${fusions}".endsWith(".junction") ? "-t STAR" : "${fusions}".endsWith(".txt") ? "-t MapSplice" : "${fusions}".endsWith(".bam") ? "-t BWA" : "-t segemehl"
    """
    CIRCexplorer2 \\
        parse \\
        $aligner \\
        $fusions \\
        -b ${prefix}.junction.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circexplorer2: \$( echo \$(CIRCexplorer2 --version 2>&1) )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junction.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circexplorer2: \$( echo \$(CIRCexplorer2 --version 2>&1) )
    END_VERSIONS
    """
}
