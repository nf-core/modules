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

process CHECKV_DOWNLOADDATABASE {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::checkv=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    input:
    path fasta
    path db

    output:
    path "checkv-db-*"           , emit: checkv_db
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def checkv_db = db ?: ''
    def update_sequence = fasta ?: ''

    // determine which command needs to be executed
    if (checkv_db != '' & update_sequence != '') {
        method = "checkv update_database --threads $task.cpus"
    }else if (checkv_db != '' & update_sequence == '') {
        method = 'ln -s'
    }else{
        method = 'checkv download_database'
    }

    """
    $method \\
        $args \\
        $db \\
        ./  \\
        $fasta \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """

}
