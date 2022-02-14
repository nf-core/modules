def VERSION = '0.3.2' // Version information not provided by tool on CLI

process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::r-ichorcna=0.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.3.2--r41hdfd78af_0' :
        'quay.io/biocontainers/r-ichorcna:0.3.2--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(wig)
    path gc_wig
    path map_wig
    path panel_of_normals
    path centromere

    output:
    tuple val(meta), path("*.cna.seg")    , emit: cna_seg
    tuple val(meta), path("*.params.txt") , emit: ichorcna_params
    path "**/*genomeWide.pdf"             , emit: genome_plot
    path "versions.yml"                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''

    """
    runIchorCNA.R --id ${prefix} \\
        $args \\
        --WIG ${wig} \\
        --id ${meta.id} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}
