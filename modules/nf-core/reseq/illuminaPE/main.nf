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
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process RESEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/reseq:1.1--py38hc35fec1_3':
        'biocontainers/reseq:1.1--py38hc35fec1_3' }"


    input:
    path(fasta)
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz")         , emit: fastqgz
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reseq illuminaPE \\
        ${args} \\
        -j ${task.cpus} \\
        -r ${fasta} \\
        -b ${bam} \\
        -1 ${prefix}.1.fastq.gz \\
        -2 ${prefix}.2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reseq: \$(reseq --version |& sed 's/ReSeq version //g')
    END_VERSIONS
    """
}
