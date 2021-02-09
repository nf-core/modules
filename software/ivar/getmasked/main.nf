// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IVAR_GETMASKED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(filtered_tsv)
    path (bed)
    path (primer_pairs)

    output:
    tuple val(meta), path("*.masked.txt"), emit: masked_sites
    path "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ivar getmasked \\
        $options.args \\
        -p $prefix \\
        -i $filtered_tsv \\
        -b $bed \\
        -f $primer_pairs

    ivar version | head -n1 2>&1 | sed 's/^.*iVar version //g' > ${software}.version.txt
    """
}
