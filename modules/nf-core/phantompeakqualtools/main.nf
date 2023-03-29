process PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::phantompeakqualtools=1.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phantompeakqualtools:1.2.2--0' :
        'quay.io/biocontainers/phantompeakqualtools:1.2.2--0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.out")  , emit: spp
    tuple val(meta), path("*.pdf")  , emit: pdf
    tuple val(meta), path("*.Rdata"), emit: rdata
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    RUN_SPP=`which run_spp.R`
    Rscript $args -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="$bam" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phantompeakqualtools: $VERSION
    END_VERSIONS
    """
}
