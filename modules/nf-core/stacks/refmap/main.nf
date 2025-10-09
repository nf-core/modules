// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process STACKS_REFMAP {
    tag '$bam'
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stacks:2.68--h077b44d_0':
        'biocontainers/stacks:2.68--h077b44d_0' }"


    input:// TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    path bam_path
    path popmap

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels


    path "catalog.calls" , emit: catalog_calls
    path "catalog.chrs.tsv" , emit: catalog_chrs
    path "catalog.fa.gz" , emit: catalog_fa
    path "gstacks.log" , emit: gstacks_log
    path "gstacks.log.distribs", emit: gstacks_log_distribs
    path "populations.haplotypes.tsv", emit: haplotypes
    path "populations.hapstats.tsv", emit: hapstats
    path "populations.sumstats.tsv", emit: sumstats
    path "populations.sumstats_summary.tsv", emit: sumstats_summary
    path "populations.log", emit: populations_log
    path "populations.log.distribs", emit: populations_log_distribs
    path "ref_map.log", emit: ref_map_log
    path "populations.snps.vcf", emit: vcf , optional: true
    path "populations.snps.genepop", emit: genepop , optional: true
    path "populations.structure", emit: structure , optional: true
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "refmap_output"
    //def bam_path = "--samples ${bam_path}" ?: ''
    //def popmap = "--popmap ${popmap}" ?: ''

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    ref_map.pl \\
        --samples ./ \\
        --popmap ${popmap} \\
        -o . \\
        -T $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stacks: \$(populations -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "refmap_output"

    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """

    touch catalog.calls
    touch catalog.chrs.tsv
    touch catalog.fa.gz
    touch gstacks.log
    touch gstacks.log.distribs
    touch populations.haplotypes.tsv
    touch populations.hapstats.tsv
    touch populations.sumstats.tsv
    touch populations.sumstats_summary.tsv
    touch populations.log
    touch populations.log.distribs
    touch ref_map.log
    touch populations.snps.vcf
    touch populations.snps.genepop
    touch populations.structure

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stacks: \$(populations -v)
    END_VERSIONS
    """
}
