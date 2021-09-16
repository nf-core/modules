// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DIAMOND_BLASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container "quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0"
    }

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.diamond_blastp.txt'), emit: txt
    path '*.version.txt', emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    DB=`find -L ./ -name "*.dmnd" | sed 's/.dmnd//'`
    diamond blastp \\
        --threads $task.cpus \\
        --db \$DB \\
        --query $fasta \\
        $options.args \\
        --out ${prefix}.diamond_blastp.txt
    echo \$(diamond --version 2>&1) | tail -n 1 | sed 's/^diamond version //' > ${software}.version.txt
    """
}
