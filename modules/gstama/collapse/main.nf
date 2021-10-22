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

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0.1--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0.1--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_tc.bed")                    , emit: bed
    tuple val(meta), path("*_tc_trans_read.bed")         , emit: bed_trans_reads
    tuple val(meta), path("*_tc_local_density_error.txt"), emit: local_density_error
    tuple val(meta), path("*_tc_polya.txt")              , emit: polya
    tuple val(meta), path("*_tc_read.txt")               , emit: read
    tuple val(meta), path("*_tc_strand_check.txt")       , emit: strand_check
    tuple val(meta), path("*_tc_trans_report.txt")       , emit: trans_report
    path "versions.yml"                                  , emit: versions

    tuple val(meta), path("*_tc_varcov.txt.txt")         , optional: true, emit: varcov
    tuple val(meta), path("*_tc_variants.txt")           , optional: true, emit: variants

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    tama_collapse.py \\
        -s $bam \\
        -f $fasta \\
        -p ${prefix}_tc \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( tama_collapse.py -version | grep 'tc_version_date_'|sed 's/tc_version_date_//g' )
    END_VERSIONS
    """
}
