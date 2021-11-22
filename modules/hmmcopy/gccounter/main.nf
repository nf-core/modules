def VERSION = '0.1.1'

process HMMCOPY_GCCOUNTER {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_5' :
        'quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_5' }"

    input:
    path fasta

    output:
    path "*.gc.wig"    , emit: wig
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    gcCounter \\
        $args \\
        ${fasta} > ${fasta.baseName}.gc.wig

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        hmmcopy: \$(echo $VERSION)
    END_VERSIONS
    """
}
