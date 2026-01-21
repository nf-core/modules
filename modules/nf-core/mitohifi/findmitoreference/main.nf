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
    tuple val(meta), path("*.fasta"), path("*.gb"), emit: reference
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    // old version command: \$(mitohifi.py -v | sed 's/.* //')
    tuple val("${task.process}"), val('mitohifi'), eval('echo 3.2.3'), emit: versions_mitohifi, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }

    def args         = task.ext.args ?: ''
    def ncbi_api_key = secrets.NCBI_API_KEY ? "--ncbi-api-key \$NCBI_API_KEY" : ""
    """
    findMitoReference.py \\
        ${ncbi_api_key} \\
        --species "$species" \\
        --outfolder . \\
        $args
    """

    stub:
    """
    touch accession.fasta
    touch accession.gb
    """
}
