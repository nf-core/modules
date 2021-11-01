// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.99.0'

process HMMCOPY_GCCOUNTER {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_5"
    } else {
        container "quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_5"
    }

    input:
    path fasta

    output:
    path "*.gc.wig"          , emit: wig
    path "versions.yml" , emit: versions

    script:
    """
    gcCounter \\
        $options.args \\
        ${fasta} > ${fasta.baseName}.gc.wig

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
