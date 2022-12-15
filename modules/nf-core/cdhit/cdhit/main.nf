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

process CDHIT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cd-hit=4.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit%3A4.8.1--h5b5514e_7':
        'quay.io/biocontainers/cd-hit:4.8.1--h5b5514e_7' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path(".txt"),    emit: clusters
    path "versions.yml"           ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threshold = params.cdhit_threshold ?: 0.9
    def wordsize = params.cdhit_wordsize ?: 5
    """
    cd-hit \\
    -i $sequences \\
    -o $prefix \\
    -c $threshold \\
    -n $wordsize \\
    -M $task.memory.bytes \\
    -T $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(echo \$(cd-hit -h|head -n 1 2>&1) | sed 's/====== CD-HIT version //;s/ (built on .*) ======//' ))
    END_VERSIONS
    """
}
