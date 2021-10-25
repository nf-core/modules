// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

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
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

process CHECKM_LINEAGEWF {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::checkm-genome=1.1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/checkm-genome:1.1.3--py_1"
    } else {
        container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"
    }

    input:
    tuple val(meta), path(fasta)
    val(fasta_ext)

    output:
    tuple val(meta), path("${prefix}"), emit: checkm
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    checkm \\
        lineage_wf \\
        -t $task.cpus \\
        --pplacer_threads $task.cpus \\
        -x $fasta_ext \\
        $options.args \\
        . \\
        ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$( checkm 2>&1 \) | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
