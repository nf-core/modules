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

process BISCUIT_INDEX {
    tag "$fasta"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.1.20211018" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.0.1.20211018--h81a5ba2_1':
        'quay.io/biocontainers/biscuit:1.0.1.20211018--h81a5ba2_1' }"

    input:
    path fasta, stageAs: "BiscuitIndex/*"

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    // TODO nf-core: List additional required output channels/values here
    path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    biscuit \\
        index \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$(echo \$(biscuit version 2>&1) | sed 's/^.*BISCUIT Version: //; s/Using.*\$//')
    END_VERSIONS
    """
}
