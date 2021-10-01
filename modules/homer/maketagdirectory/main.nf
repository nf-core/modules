// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '4.11'

process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3"
    } else {
        container "quay.io/biocontainers/homer:4.11--pl526hc9558a2_3"
    }

    input:
    tuple val(meta), path(bed)
    path fasta

    output:
    tuple val(meta), path("tag_dir"), emit: tagdir
    path  "versions.yml"             , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    makeTagDirectory \\
        tag_dir \\
        $options.args \\
        $bed \\
        -genome $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
