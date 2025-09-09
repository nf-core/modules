process METAMAPS_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metamaps:0.1.633d2e0--h21ec9f0_0':
        'biocontainers/metamaps:0.1.633d2e0--h21ec9f0_0' }"

    input:
    tuple val(meta), path(classification_res), path(meta_file), path(meta_unmappedreadsLengths), path(para_file)
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
        $args \\
        --mappings $classification_res \\
        --threads $task.cpus \\
        --DB $database_folder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps: \$(metamaps | sed -n 2p | sed 's/^.*MetaMaps v //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_classification_res.EM.WIMP
    touch ${prefix}_classification_res.EM.evidenceUnknownSpecies
    touch ${prefix}_classification_res.EM.reads2Taxon
    touch ${prefix}_classification_res.EM
    touch ${prefix}_classification_res.EM.contigCoverage
    touch ${prefix}_classification_res.EM.lengthAndIdentitiesPerMappingUnit
    touch ${prefix}_classification_res.EM.reads2Taxon.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps: \$(metamaps | sed -n 2p | sed 's/^.*MetaMaps v //')
    END_VERSIONS
    """

}
