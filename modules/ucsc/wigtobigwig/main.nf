// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '377' // No version information printed

process UCSC_WIGTOBIGWIG {
    tag '$wig'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::ucsc-wigtobigwig=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-wigtobigwig:377--h0b8a92a_2"
    } else {
        container "quay.io/biocontainers/ucsc-wigtobigwig:377--h0b8a92a_2"
    }

    input:
    path wig
    path chromsizes

    output:
    path "*.bw"                   , emit: bw
    path "versions.yml"           , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    wigToBigWig \\
        $options.args \\
        $wig \\
        $chromsizes \\
        ${wig.getSimpleName()}.bw

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo "$VERSION")
    END_VERSIONS
    """
}
