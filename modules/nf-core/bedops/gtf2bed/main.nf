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

process BEDOPS_GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0':
        'biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    input:
    path gtf

    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    cat \\
    $gtf \\
    | gtf2bed \\
    $args \\
    --attribute-key=exon_id \\
    > {prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    bedops: \$(bedops --version | grep version | awk ' { print $2 } ')
    END_VERSIONS
    """
}
