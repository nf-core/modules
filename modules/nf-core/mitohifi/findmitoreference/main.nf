process MITOHIFI_FINDMITOREFERENCE {
    tag "$species"
    label 'process_single'
    secret secrets.NCBI_API_KEY ? "NCBI_API_KEY" : ""

    // NOTE: An optional NCBI API key can be supplied to MITOHIFI_FINDMITOREFERENCE.
    // This should be set using Nextflow's secrets functionality:
    // `nextflow secrets set NCBI_API_KEY <key>`
    //
    // See https://www.nextflow.io/docs/latest/secrets.html for more information.

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }

    def args         = task.ext.args ?: ''
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION      = '3.2.3'
    def ncbi_api_key = secrets.NCBI_API_KEY ? "--ncbi-api-key \$NCBI_API_KEY" : ""
    """
    findMitoReference.py \\
        ${ncbi_api_key} \\
        --species "$species" \\
        --outfolder . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: ${VERSION}
    END_VERSIONS
    """

    stub:
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '3.2.3'
    """
    touch accession.fasta
    touch accession.gb

    ## old version command: \$(mitohifi.py -v | sed 's/.* //')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: ${VERSION}
    END_VERSIONS
    """
}
