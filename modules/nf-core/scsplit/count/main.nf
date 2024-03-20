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

process SCSPLIT_COUNT {
    tag "$meta.id"
    label 'process_high'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/irenerobles93/scsplit:1.0.8"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //    'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(barcode)

    output:
    tuple val(meta), path("*alt_filtered.csv"), emit: alf_filtered
    tuple val(meta), path("*ref_filtered.csv"), emit: ref_filtered
    tuple val(meta), path("*scSplit.log")     , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scSplit count -v $vcf -i $bam -b $barcode -r ${prefix}_ref_filtered.csv -a ${prefix}_alt_filtered.csv -o .
    mv scSplit.log ${prefix}_scSplit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scsplit: 1.0.8
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ref_filtered.csv
    touch ${prefix}_alt_filtered.csv
    touch ${prefix}_scSplit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scsplit: 1.0.8
    END_VERSIONS
    """
}
