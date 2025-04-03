process WISECONDORX_PREDICT {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(npz)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(blacklist)

    output:
    tuple val(meta), path("*_aberrations.bed")      , emit: aberrations_bed, optional:true
    tuple val(meta), path("*_bins.bed")             , emit: bins_bed, optional:true
    tuple val(meta), path("*_segments.bed")         , emit: segments_bed, optional:true
    tuple val(meta), path("*_statistics.txt")       , emit: chr_statistics, optional:true
    tuple val(meta), path("[!genome_wide]*.png")    , emit: chr_plots, optional:true
    tuple val(meta), path("genome_wide.png")        , emit: genome_plot, optional:true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--bed --plot'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = blacklist ? "--blacklist ${blacklist}" : ""

    def plots = args.contains("--plot") ? "mv ${prefix}.plots/* ." : ""

    def VERSION = '1.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    WisecondorX predict \\
        ${npz} \\
        ${reference} \\
        ${prefix} \\
        ${bed} \\
        ${args}

    ${plots}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '--bed --plot'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def bed = args.contains("--bed") ? "touch ${prefix}_aberrations.bed && touch ${prefix}_bins.bed && touch ${prefix}_statistics.txt && touch ${prefix}_segments.bed" : ""
    def plot = args.contains("--plot") ? "touch genome_wide.png && touch chr22.png && touch chr1.png" : ""

    """
    ${bed}
    ${plot}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """
}
