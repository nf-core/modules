#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
def SMNCOPYNUMBERCALLER_VERSION = 'SMNCopyNumberCaller commit 3e67e3b on Feb 8, 2020'
// No versioning included with this program; added github commit from:
// https://github.com/Illumina/SMNCopyNumberCaller

process SMNCOPYNUMBERCALLER {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used with SMNCopyNumberCaller at the moment. Please use Docker or Singularity containers."
    }
    container "clinicalgenomics/smncopynumbercaller:v1.1.2"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("manifest.txt"), emit: manifest
    tuple val(meta), path("*.tsv"), emit: smncopynumber
    tuple val(meta), path("*.json"), emit: run_metrics
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    def manifest_in = file(task.workDir + 'manifest.txt')
    manifest_in.text = bam.join("\n")

    script:
    def args = task.ext.args ?: ''
    def cpus = task.cpus
    def genome_version = task.ext.genome_version // [19/37/38]
    def out_dir = task.outDir
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    smn_caller.py \\
        $args \\
        --manifest $manifest_in \\
        --genome $genome_version \\
        --prefix $prefix \\
        --outDir $out_dir \\
        --threads $cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SMNCopyNumberCaller: $SMNCOPYNUMBERCALLER_VERSION
    END_VERSIONS
    """
}
