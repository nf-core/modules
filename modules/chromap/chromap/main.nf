// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

def VERSION = 0.1 // No version information printed

process CHROMAP_CHROMAP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::chromap=0.1 bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0"
    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(reads)
    path fasta
    path index
    path barcodes
    path whitelist
    path chr_order
    path pairs_chr_order

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.bed.gz")     , optional:true, emit: bed
    tuple val(meta), path("*.bam")        , optional:true, emit: bam
    tuple val(meta), path("*.tagAlign.gz"), optional:true, emit: tagAlign
    tuple val(meta), path("*.pairs.gz")   , optional:true, emit: pairs
    // TODO nf-core: List additional required output channels/values here
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def args     = options.args.tokenize()

    def file_extension = options.args.contains("--SAM")? 'sam' :
                        options.args.contains("--TagAlign")? 'tagAlign' :
                        options.args.contains("--pairs")? 'pairs' : 'bed'
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$options.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
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
            -t ${task.cpus} \\
            -x $index \\
            -r $fasta \\
            -1 ${reads.join(',')} \\
            -o ${prefix}.${file_extension}
        echo "$VERSION" > ${software}.version.txt
        """ + compression_cmds
    } else {
        """
        chromap ${args.join(' ')} \\
            -t ${task.cpus} \\
            -x $index \\
            -r $fasta \\
            -1 reads[0] \\
            -2 reads[1] \\
            -o ${prefix}.${file_extension}
        echo "$VERSION" > ${software}.version.txt
        """ + compression_cmds
    }
}
