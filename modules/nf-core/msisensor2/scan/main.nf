process MSISENSOR2_SCAN {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor2:0.1--hd03093a_0':
        'biocontainers/msisensor2:0.1--hd03093a_0' }"

    input:
    path fasta
    val output

    output:
    path output_path   , emit: scan
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def inputs      = fasta.collect{ "-d $it"}.join(" ")
    output_path = output ?: "output.scan"
    """
    msisensor2 scan \\
        $args \\
        $inputs \\
        -o $output_path

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """

    stub:
    output_path = output ?: "output.scan"
    """
    touch $output_path

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """
}
