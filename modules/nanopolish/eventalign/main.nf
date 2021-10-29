// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NANOPOLISH_EVENTALIGN {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::nanopolish==0.13.2" : null)
    container "quay.io/biocontainers/nanopolish:0.13.2--he3b7ca5_2"

    input:
    tuple val(meta), path(bam), path(bai)
    path fast5
    path fastq
    path fasta
    path gtf

    output:
    tuple val(meta), path("*eventalign.txt"), path("*summary.txt"), emit: nanopolish_outputs
    path "*.version.txt"     ,emit: version

    script:
    sample_summary = "$meta.id" +"_summary.txt"
    sample_eventalign = "$meta.id" +"_eventalign.txt"
    """
    nanopolish index -d $fast5 $fastq
    nanopolish eventalign  --reads $fastq --bam $bam --genome $fasta --scale-events --signal-index --summary $sample_summary --threads $threads > $sample_eventalign
    nanopolish --version | sed -e "s/nanopolish version //g" | head -n 1 > nanopolish.version.txt
    """
}
