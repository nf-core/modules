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

process ARCASHLA_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::arcas-hla=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.5.0--hdfd78af_1':
        'biocontainers/arcas-hla:0.5.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(fq)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.json"), emit: genotype
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.5.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    arcasHLA reference --rebuild
    arcasHLA reference --version

    arcasHLA \\
        genotype \\
        $args \\
        -t $task.cpus \\
        -o ${prefix} \\
        $fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: $VERSION
    END_VERSIONS
    """
}
