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

process FREYJA_BOOT {
    tag "$meta.id"
    label 'process_long'

    conda "bioconda::freyja=1.3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.3.12--pyhdfd78af_0':
        'quay.io/biocontainers/freyja:1.3.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(variants)
    tuple val(meta2), path(depths)
    val repeats
    path barcodes
    path lineages

    output:
    tuple val(meta), path("*.lineages.csv")     , emit: lineages
    tuple val(meta), path("*.summarized.csv")   , emit: summarized
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def repeats= repeats.isInteger() ? "--nb $repeats" : ""

    """
    freyja \\
        boot \\
        $args \\
        --nt $task.cpus \\
        $repeats \\
        --output_base $prefix \\
        --barcodes $barcodes \\
        --meta $lineages \\
        $variants \\
        $depths

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(freyja --version 2>&1) | sed 's/^.*version//' ))
    END_VERSIONS
    """
}
