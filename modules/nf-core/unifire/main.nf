process UNIFIRE {
    tag "$meta.id"
    label 'process_medium'

    container "dockerhub.ebi.ac.uk/uniprot-public/unifire:2025.1" // TODO: Update once Bioconda is available
    containerOptions {
        if (workflow.containerEngine in ['singularity', 'apptainer']) {
            return "--bind unifire:/volume"
        } else {
            return "-v ./unifire:/volume"
        }
    }

    input:
    tuple val(meta), path(faa, stageAs: "unifire/proteins.fasta")

    output:
    tuple val(meta), path("unifire/predictions_arba.out")         , emit: arba
    tuple val(meta), path("unifire/predictions_unirule.out")      , emit: unirule
    tuple val(meta), path("unifire/predictions_unirule-pirsr.out"), emit: pirsr
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2025.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # This tool needs a specific folder to be mounted to work.
    # Run UniFIRE workflow
    /opt/scripts/bin/unifire-workflow.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UniFIRE: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2025.1'
    """
    mkdir -p unifire
    touch unifire/predictions_arba.out
    touch unifire/predictions_unirule.out
    touch unifire/predictions_unirule-pirsr.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UniFIRE: ${VERSION}
    END_VERSIONS
    """
}
