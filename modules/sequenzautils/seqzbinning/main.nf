// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQUENZAUTILS_SEQZBINNING {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h8f06ca0_4"
    } else {
        container "quay.io/biocontainers/sequenza-utils:3.0.0--py27h304d29a_4"
    }

    input:
    tuple val(meta), path(seqz)

    output:
    tuple val(meta), path("*.binned.gz"), emit: seqz_binned
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    sequenza-utils \\
        seqz_binning \\
        --seqz ${seqz} \\
        $options.args \\
        -o ${prefix}.binned.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
