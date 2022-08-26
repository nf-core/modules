process MSISENSOR2_SCAN {
    tag '$fasta'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::msisensor2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/msisensor2:0.1_cv1' }"

    input:
    path fasta
    val output

    output:
    path output, emit: scan
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def inputs      = fasta.collect{ "-d $it"}.join(" ")
    """
    msisensor2 scan \\
        $args \\
        $inputs \\
        -o $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """
}
