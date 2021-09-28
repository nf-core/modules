// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PBBAM_PBMERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pbbam=1.7.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbbam:1.7.0--h058f120_1"
    } else {
        container "quay.io/biocontainers/pbbam:1.7.0--h058f120_1"
    }

    input:
    tuple val(meta), path("*.bam")

    output:
    tuple val(meta), path("*.m.bam")  , emit: bam
    tuple val(meta), path("*.bam.pbi"), emit: pbi
    path "versions.yml"               , emit: version

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pbmerge -o ${prefix}.m.bam $options.args *.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        pbbam/pbmerge: \$( pbmerge --version|sed 's/pbmerge //' )
    END_VERSIONS
    """
}
