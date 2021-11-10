// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PARACLU {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::paraclu=10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/paraclu%3A10--h9a82719_1"
    } else {
        container "quay.io/biocontainers/paraclu:10--h9a82719_1"
    }

    input:
    tuple val(meta), path(bed)
    val(min_cluster)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def VERSION=10
    """

    awk -F "\t" '{print\$1"\t"\$6"\t"\$2"\t"\$5}' < $bed > ${bed}_4P
    sort -k1,1 -k3n ${bed}_4P > ${bed}_4Ps
    paraclu $min_cluster ${bed}_4Ps > ${prefix}.clustered
    paraclu-cut  ${prefix}.clustered >  ${prefix}.clustered.simplified
    awk -F '\t' '{print \$1"\t"\$3"\t"\$4"\t"\$1":"\$3".."\$4","\$2"\t"\$6"\t"\$2}' ${prefix}.clustered.simplified >  ${prefix}.clustered.simplified.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: $VERSION
    END_VERSIONS
    """
}
