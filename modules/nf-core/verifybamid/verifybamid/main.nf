// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules
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

process VERIFYBAMID_VERIFYBAMID {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::verifybamid=1.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid%3A1.1.3--h5b5514e_6':
        'quay.io/biocontainers/verifybamid:1.1.3--h5b5514e_6' }"

    input:
    tuple val(meta), path(bam)
    path refvcf

    output:
    tuple val(meta), path("*.log")                   , optional:true, emit: log
    tuple val(meta), path("*.selfSM")                , optional:true, emit: selfsm
    tuple val(meta), path("*.depthSM")               , optional:true, emit: depthsm
    tuple val(meta), path("*.selfRG")                , optional:true, emit: selfrg
    tuple val(meta), path("*.depthRG")               , optional:true, emit: depthrg
    tuple val(meta), path("*.bestSM")                , optional:true, emit: bestsm
    tuple val(meta), path("*.bestRG")                , optional:true, emit: bestrg
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()

    def bam_file = "${bam}.endsWith('.bam')" ? "--bam ${bam}" : ""
    def refvcf_args = "${refvcf}".endsWith(".vcf") ? "--vcf ${refvcf}" : ""
    def out_args = "--out ${prefix}"

    """
    verifyBamID \\
        ${bam_file} \\
        ${refvcf_args} \\
        ${out_args} \\
        ${args_list.join(' ')} \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid: \$(echo \$(verifyBamID --help 2>&1 | sed 's/^verifyBamID //'))
    END_VERSIONS
    """
}
