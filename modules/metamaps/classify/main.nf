process METAMAPS_CLASSIFY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::metamaps=0.1.98102e9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metamaps:0.1.98102e9--h176a8bc_0':
        'quay.io/biocontainers/metamaps:0.1.98102e9--h176a8bc_0' }"

    input:
    tuple val(meta), path(classification_res)
    path database_folder

    output:
    tuple val(meta), path("*classification_res.EM.WIMP")                              , emit: wimp
    tuple val(meta), path("*classification_res.EM.evidenceUnknownSpecies")            , emit: evidence_unknown_species
    tuple val(meta), path("*classification_res.EM.reads2Taxon")                       , emit: reads2taxon
    tuple val(meta), path("*classification_res.EM")                                   , emit: em
    tuple val(meta), path("*classification_res.EM.contigCoverage")                    , emit: contig_coverage
    tuple val(meta), path("*classification_res.EM.lengthAndIdentitiesPerMappingUnit") , emit: length_and_id
    tuple val(meta), path("*classification_res.EM.reads2Taxon.krona")                 , emit: krona
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metamaps \\
        classify \\
        $args
        --mappings $classification_res \\
        --threads $task.cpus \\
        --DB $database_folder \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps: \$(metamaps | sed -n 2p | sed 's/^.*MetaMaps v //')
    END_VERSIONS
    """
}
