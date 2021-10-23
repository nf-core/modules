// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA_MERGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0.2--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0.2--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bed)
    path filelist

    output:
    tuple val(meta), path("*.bed")             , emit: bed
    tuple val(meta), path("*_gene_report.txt") , emit: gene_report
    tuple val(meta), path("*_merge.txt")       , emit: merge
    tuple val(meta), path("*_trans_report.txt"), emit: trans_report
    path "versions.yml"                        , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    tama_merge.py \\
        -f $filelist \\
        -d merge_dup \\
        -p ${prefix} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( tama_merge.py -version | head -n1 )
    END_VERSIONS
    """
}
