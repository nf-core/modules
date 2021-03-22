// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/software
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join

// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, "fake files" MAY be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::adapterremoval=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/adapterremoval:2.3.2--hb7ba0dd_0"
    } else {
        container "quay.io/biocontainers/adapterremoval:2.3.2--hb7ba0dd_0"
    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(reads)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log'), emit: log
    // TODO nf-core: List additional required output channels/values here
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (meta.single_end) {
        """
        AdapterRemoval  \\
            --file1 $reads \\
            $options.args \\
            --basename $prefix \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --output1 ${prefix}.trimmed.fastq.gz \\
            --gzip \\

        AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g" > ${software}.version.txt
        """
    } else if (!meta.single_end && !meta.collapse) {
        """
        AdapterRemoval  \\
            --file1 ${reads[0]} \\
            --file2 ${reads[0]} \\
            $options.args \\
            --basename $prefix \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --output1 ${prefix}.pair1.trimmed.fastq.gz \\
            --output2 ${prefix}.pair2.trimmed.fastq.gz \\
            --gzip \\

        AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g" > ${software}.version.txt
        """
    } else {
        """
        AdapterRemoval  \\
            --file1 ${reads[0]} \\
            --file2 ${reads[0]} \\
            --collapse \\
            $options.args \\
            --basename $prefix \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --gzip \\
        
        cat *.collapsed.gz *.collapsed.truncated.gz > ${prefix}.merged.fastq.gz
        AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g" > ${software}.version.txt
        """
    }
}
