process PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a1cddfad5b503ee347cc5de17d172e1876c547fca00aa844559c9e764fb400f/data' :
        'community.wave.seqera.io/library/phantompeakqualtools:1.2.2--f8026fe2526a5e18' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.out")  , emit: spp
    tuple val(meta), path("*.pdf")  , emit: pdf
    tuple val(meta), path("*.Rdata"), emit: rdata
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), topic: versions, emit: versions_r
    tuple val("${task.process}"), val('phantompeakqualtools'), val('1.2.2'), topic: versions, emit: versions_phantompeakqualtools // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    RUN_SPP=`which run_spp.R`
    Rscript $args -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="$bam" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" $args2
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.spp.pdf
    touch ${prefix}.spp.Rdata
    touch ${prefix}.spp.out
    """
}
