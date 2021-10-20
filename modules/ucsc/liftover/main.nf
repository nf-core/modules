// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = 'v377'

process UCSC_LIFTOVER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ucsc-liftover=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-liftover:377--h0b8a92a_3"
    } else {
        container "quay.io/biocontainers/ucsc-liftover:377--h0b8a92a_3"
    }

    input:
    tuple val(meta), path(bed)
    path(chain)

    output:
    tuple val(meta), path("*.lifted.bed")  , emit: lifted
    tuple val(meta), path("*.unlifted.bed"), emit: unlifted
    path "versions.yml"                    , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    liftOver \\
        $options.args \
        $bed \\
        $chain \\
        ${prefix}.lifted.bed \\
        ${prefix}.unlifted.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo "$VERSION")
    END_VERSIONS
    """
}
