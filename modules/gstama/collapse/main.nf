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
    path genome

    output:
    tuple val(meta), path("*_tc.bed")                    , emit: bed
    tuple val(meta), path("*_tc_trans_read.bed")         , emit: bed_trans_reads
    tuple val(meta), path("*_tc_local_density_error.txt"), emit: local_density_error
    tuple val(meta), path("*_tc_polya.txt")              , emit: tc_polya
    tuple val(meta), path("*_tc_read.txt")               , emit: tc_read
    tuple val(meta), path("*_tc_strand_check.txt")       , emit: tc_strand_check
    tuple val(meta), path("*_tc_trans_report.txt")       , emit: tc_trans_report
    path "versions.yml"                                  , emit: version

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    NLINES=\$(samtools view $bam|wc -l)
    MODE=""

    if [ "\$NLINES" -gt 1000 ]; then
        MODE="-rm low_mem"
    fi

    tama_collapse.py \\
        -s $bam \\
        -f $genome \\
        -p ${prefix}_tc \\
        \$MODE \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        tama collapse: \$( tama_collapse.py -version | grep 'tc_version_date_'|sed 's/tc_version_date_//g' )
    END_VERSIONS
    """
}
