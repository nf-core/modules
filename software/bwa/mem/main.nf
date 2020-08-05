// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    //container "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"

    conda (params.conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(reads)
    path index
    path fasta
    val options

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def rg       = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    bwa mem \\
        $ioptions.args \\
        $rg \\
        -t $task.cpus \\
        $fasta \\
        $reads \\
        | samtools view $ioptions.args2 -@ $task.cpus -bS -o ${prefix}.bam -

    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
}
