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


process MACSE_REFINEALIGNMENT {
    tag "$meta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macse:2.07--hdfd78af_0':
        'biocontainers/macse:2.07--hdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("*_NT"), emit: result_fas_NT
    tuple val(meta), path("*_AA"), emit: result_fas_AA
    tuple val("${task.process}"), val("macse"), eval("macse --version 2>&1 | sed 's/ (.*) //g'"), topic: versions, emit: versions_macse

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName

    """
    macse -prog refineAlignment \
        -align ${fasta} \
        $args \\
        -out_NT ${prefix}_NT \\
        -out_AA ${prefix}_AA \\
    """

    stub:
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    echo $args
    touch ${prefix}.fas
    """
}
