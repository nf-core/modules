// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MLST {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::mlst=2.19.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mlst:2.19.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/mlst:2.19.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mlst \\
      --threads $task.cpus \\
      $fasta \\
      > ${fasta}.mlst.tsv

    mlst --version > ${software}.version.txt

    echo "MLST" > versions.yml
    echo "  mlst: \$(mlst --version | grep -o '[0-9.]*')" >> versions.yml

    """
    //echo "${getProcessName(task.process)}:" > versions.yml

    //samtools \\
    //    sort \\
    //    $options.args \\
    //    -@ $task.cpus \\
    //    -o ${prefix}.bam \\
    //    -T $prefix \\
    //    $bam

    //echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt

}
