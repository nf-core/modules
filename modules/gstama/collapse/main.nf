// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA_COLLAPSE {
    tag "$meta.id"
    label 'process_tama'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)
    path genome

    output:
    tuple val(meta), path("*_tc.bed")               , emit: bed
    tuple path("*_tc_trans_read.bed"), path("*.txt"), emit: reports
    path "*.version.txt"                            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = bam.toString().replaceAll(/.bam/, "") + "_tc"
    """
    nlines=\$(samtools view $bam|wc -l)
    mode=""

    if [ \$nlines -gt 1000 ]; then
	    mode="-rm low_mem"
    fi

    tama_collapse.py \\
        -s $bam \\
        -f $genome \\
        -p $prefix \\
        \$mode \\
        $options.args

    echo \$(tama_collapse.py -version 2>&1) | grep 'tc_version_date_' > ${software}.version.txt
    """
}
