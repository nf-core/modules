process MSISENSOR2_MSI {
    tag "$meta.id"
    label 'process_med'

    conda (params.enable_conda ? "bioconda::msisensor2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/msisensor2:0.1_cv1' }"

    input:
    tuple val(meta), path(tumour_bam), path(normal_bam), path(intervals)
    path models
    path

    output:
    tuple val(meta), path("${prefix}")        , emit: msi
    tuple val(meta), path("${prefix}_dis")    , emit: distribution
    tuple val(meta), path("${prefix}_somatic"), emit: somatic
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumour_bam_cmd   = tumour_bam ? "-t $tumour_bam" : ""
    def normal_bam_cmd   = normal_bam ? "-n $normal_bam" : ""
    def interval_command = intervals  ? "-e $intervals"  : ""
    """
    msisensor2 msi \\
        $args \\
        -M $models \\
        $interval_command \\
        $tumour_bam_cmd \\
        $normal_bam_cmd \\
        -o $prefix


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """
}
