process MITOHIFI_FINDMITOREFERENCE {
    tag "$species"
    label 'process_single'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    tuple val(meta), val(species)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.gb")   , emit: gb
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '3.2.3' // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    """
    findMitoReference.py \\
        --species "$species" \\
        --outfolder . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: ${VERSION}
    END_VERSIONS
    """

    stub:
    """
    touch accession.fasta
    touch accession.gb

    ## old version command: \$(mitohifi.py -v | sed 's/.* //')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
    """
}
