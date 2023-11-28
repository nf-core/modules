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

process CHEWBBACA_ALLELECALL {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::chewbbaca=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.3.1--pyhdfd78af_0' :
        'biocontainers/chewbbaca:3.3.1--pyhdfh78af_0' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(fasta)
    path(scheme)

    output:

    tuple val(meta), path("*_results_statistics.tsv"),      emit: stats
    tuple val(meta), path("*_results_contigsInfo.tsv"),     emit: contigsInfo
    tuple val(meta), path("*_results_alleles.tsv"),         emit: alleles
    tuple val(meta), path("*_paralogous_counts.tsv"),       emit: paralogous_counts, optional:true
    tuple val(meta), path("*_paralogous_loci.tsv"),         emit: paralogous_loci, optional:true
    tuple val(meta), path("*_logging_info.txt"),            emit: log
    tuple val(meta), path("*_cds_coordinates.tsv"),         emit: cds_coordinates, optional:true
    tuple val(meta), path("*_invalid_cds.txt"),             emit: invalid_cds, optional:true
    tuple val(meta), path("*_loci_summary_stats.tsv"),      emit: loci_summary_stats,   optional:true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    chewie \
      AlleleCall \
      --cpu ${task.cpus} \
      --input-files ${fasta} \
      --schema-directory ${scheme} \
      --output-directory results

    mv results/${prefix}_results_statistics.tsv ${prefix}_results_statistics.tsv
    mv results/${prefix}_results_contigsInfo.tsv ${prefix}_results_contigsInfo.tsv
    mv results/${prefix}_results_alleles.tsv ${prefix}_results_alleles.tsv
    mv results/${prefix}_paralogous_counts.tsv ${prefix}_paralogous_counts.tsv
    mv results/${prefix}_paralogous_loci.tsv ${prefix}_paralogous_loci.tsv
    mv results/${prefix}_logging_info.txt ${prefix}_logging_info.txt
    mv results/${prefix}_cds_coordinates.tsv ${prefix}_cds_coordinates.tsv
    mv results/${prefix}_invalid_cds.txt ${prefix}_invalid_cds.txt
    mv results/${prefix}_loci_summary_stats.tsv ${prefix}_loci_summary_stats.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(chewie --version 2>&1) | sed 's/^.*chewie //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}_results_statistics.tsv
    touch ${prefix}_results_contigsInfo.tsv
    touch ${prefix}_results_alleles.tsv
    touch ${prefix}_paralogous_counts.tsv
    touch ${prefix}_paralogous_loci.tsv
    touch ${prefix}_logging_info.txt
    touch ${prefix}_cds_coordinates.tsv
    touch ${prefix}_invalid_cds.txt
    touch ${prefix}_loci_summary_stats.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(chewie --version 2>&1) | sed 's/^.*chewie //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
