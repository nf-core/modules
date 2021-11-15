// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAMUTIL_TRIMBAM {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bamutil=1.0.15" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bamutil:1.0.15--h2e03b76_1"
    } else {
        container "quay.io/biocontainers/bamutil:1.0.15--h2e03b76_1"
    }

    input:
    tuple val(meta), path(bam), val(trim_left), val(trim_right)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bam \\
        trimBam \\
        $bam \\
        ${prefix}.bam \\
        $options.args \\
        -L $trim_left \\
        -R $trim_right

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$( bam trimBam 2>&1 ) | sed 's/^Version: //;s/;.*//' )
    END_VERSIONS
    """
}
