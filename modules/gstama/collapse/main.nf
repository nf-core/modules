// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA_COLLAPSE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0"

    }

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_collapsed.bed")          , emit: bed
    tuple val(meta), path("*_trans_read.bed")         , emit: bed_trans_reads
    tuple val(meta), path("*_local_density_error.txt"), emit: local_density_error
    tuple val(meta), path("*_polya.txt")              , emit: polya
    tuple val(meta), path("*_read.txt")               , emit: read
    tuple val(meta), path("*_strand_check.txt")       , emit: strand_check
    tuple val(meta), path("*_trans_report.txt")       , emit: trans_report
    path "versions.yml"                               , emit: versions

    tuple val(meta), path("*_varcov.txt")             , emit: varcov  , optional: true
    tuple val(meta), path("*_variants.txt")           , emit: variants, optional: true

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    tama_collapse.py \\
        -s $bam \\
        -f $fasta \\
        -p ${prefix} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( tama_collapse.py -version | grep 'tc_version_date_'|sed 's/tc_version_date_//g' )
    END_VERSIONS
    """
}
