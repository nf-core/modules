process COOLER_DIGEST {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    input:
    path fasta
    path chromsizes
    val  enzyme

    output:
    path "*.bed", emit: bed
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions_cooler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cooler digest \\
        $args \\
        -o "${fasta.baseName}_${enzyme.replaceAll(/[^0-9a-zA-Z]+/, '_')}.bed" \\
        $chromsizes \\
        $fasta \\
        $enzyme

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
