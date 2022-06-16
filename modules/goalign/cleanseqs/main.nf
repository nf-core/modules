process GOALIGN_CLEANSEQS {
    tag "goalign-cleanseqs:${aln}";
    label 'process_low'

    conda (params.enable_conda ? "bioconda::goalign=0.3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goalign:0.3.5--h65a6115_0':
        'quay.io/biocontainers/goalign:0.3.5--h65a6115_0' }"

    input:
    path(aln)

    output:
    path("*.{fasta,fas,fa,phy}"), emit: aln
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${aln.getName()}_clean"
    """
    goalign \\
        clean \\
        seqs \\
        -t ${task.cpus} \\
        $args \\
        -o ${prefix}.${aln.getExtension()} \\
        -i $aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goalign: \$(echo \$(goalign version))
    END_VERSIONS
    """
}
