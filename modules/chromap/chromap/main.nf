// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = 0.1 // No version information printed

process CHROMAP_CHROMAP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::chromap=0.1 bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0"
    }

    input:
    tuple val(meta), path(reads)
    path fasta
    path index
    path barcodes
    path whitelist
    path chr_order
    path pairs_chr_order

    output:
    tuple val(meta), path("*.bed.gz")     , optional:true, emit: bed
    tuple val(meta), path("*.bam")        , optional:true, emit: bam
    tuple val(meta), path("*.tagAlign.gz"), optional:true, emit: tagAlign
    tuple val(meta), path("*.pairs.gz")   , optional:true, emit: pairs
    path "*.version.txt"                                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def args     = options.args.tokenize()

    def file_extension = options.args.contains("--SAM")? 'sam' :
                        options.args.contains("--TagAlign")? 'tagAlign' :
                        options.args.contains("--pairs")? 'pairs' : 'bed'
    if (barcodes) {
        args << "-b ${barcodes.join(',')}"
        if (whitelist) {
            args << "--barcode-whitelist $whitelist"
        }
    }
    if (chr_order) {
        args << "--chr-order $chr_order"
    }
    if (pairs_chr_order){
        args << "--pairs-natural-chr-order $pairs_chr_order"
    }
    def compression_cmds = """
    gzip ${prefix}.${file_extension}
    """
    if (options.args.contains("--SAM")) {
        compression_cmds = """
        samtools view $options.args2 -@ ${task.cpus} -bh \\
            -o ${prefix}.bam ${prefix}.${file_extension}
        rm ${prefix}.${file_extension}

        samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
        """
    }
    if (meta.single_end) {
        """
        chromap ${args.join(' ')} \\
            -t $task.cpus \\
            -x $index \\
            -r $fasta \\
            -1 ${reads.join(',')} \\
            -o ${prefix}.${file_extension}
        echo "$VERSION" > ${software}.version.txt
        """ + compression_cmds
    } else {
        """
        chromap ${args.join(' ')} \\
            -t $task.cpus \\
            -x $index \\
            -r $fasta \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${prefix}.${file_extension}
        echo "$VERSION" > ${software}.version.txt
        """ + compression_cmds
    }
}
