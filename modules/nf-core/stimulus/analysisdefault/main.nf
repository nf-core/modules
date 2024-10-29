// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// DONE nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// DONE nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// DONE nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process STIMULUS_ANALYSISDEFAULT {
    tag "$meta.id"
    label 'process_medium'

    // DONE nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // DONE nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container "docker.io/mathysgrapotte/stimulus-py:latest"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    path(data)
    path(experiment_config)
    path(model_config)
    path(weights)
    path(optimizer)
    path(metrics)
    path(model)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*_modelcheck.log"), path("performance_tune_train/"), path("performance_robustness/"), emit: analysis
    // DONE nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // DONE nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // DONE nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // DONE nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // DONE nf-core: Please replace the example samtools command below with your module's command
    // DONE nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    stimulus-analysis-default \\
        $args \\
        -m ${model} \\
        -w ${weights} \\
        -me ${metrics} \\
        -ec ${experiment_config} \\
        -mc ${model_config} \\
        -d ${data} \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stimulus: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def STIMULUS_VER = '0.0.9' // container not used in stub, change manually
    """
    touch ${meta.id}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: ${STIMULUS_VER}
    END_VERSIONS
    """
}
