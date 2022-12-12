#!/usr/bin/env nextflow

process SMNCOPYNUMBERCALLER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::smncopynumbercaller=1.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smncopynumbercaller:1.1.2--py310h7cba7a3_0' :
        'quay.io/biocontainers/smncopynumbercaller:1.1.2--py310h7cba7a3_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("out/*.tsv"), emit: smncopynumber
    tuple val(meta), path("out/*.json"), emit: run_metrics
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    manifest_text = bam.join("\n")
    def args = task.ext.args ?: ''
    def genome_version = task.ext.genome_version // [19/37/38]
    def prefix = task.ext.prefix ?: "${meta.id}"
    def smncopynumbercaller_version = "1.1.2"
    """
    echo "$manifest_text" >manifest.txt
    smn_caller.py \\
        $args \\
        --manifest manifest.txt \\
        --genome $genome_version \\
        --prefix $prefix \\
        --outDir "out" \\
        --threads "${task.cpus}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SMNCopyNumberCaller: $smncopynumbercaller_version
    END_VERSIONS
    """
}
