process GOATOOLS_FINDENRICHMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2d9b75ea2491891521edf480d33b172f2c9ba68b2abae6e0b146ecf0f1d4d3bf/data':
        'community.wave.seqera.io/library/goatools_statsmodels:a6da305ad60b694a' }"

    input:
    tuple val(meta), path(study_list), path(population_list), path(annotation)
    path(obo)
    path(goslim)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('goatools'), eval("goatools --version"), topic: versions, emit: versions_goatools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def obo_arg = obo ? "--obo $obo" : ""
    def goslim_arg = goslim ? "--goslim $goslim" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    goatools \\
        find_enrichment
        $args \\
        $obo_arg \\
        $goslim_arg \\
        $study_list \\
        $population_list \\
        $annotation \\
        --outfile ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tsv
    """
}
