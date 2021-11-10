// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.1.1'

process HMMCOPY_GENERATEMAP {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // using a container generated from the hmmcopyutils github, but with a working copy of generateMap.pl, including the bowtie dependency and fastaToRead script.
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://crukmi/hmmcopyutils"
    } else {
        container "docker://crukmi/hmmcopyutils"
    }

    input:
    path fasta

    output:
    path "*.map.bw"              , emit: bigwig
    path "versions.yml"          , emit: versions

    script:

    """
    # build required indexes
    generateMap.pl -b \\
        $options.args \\
        $fasta

    # run
    generateMap.pl \\
        $options.args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
