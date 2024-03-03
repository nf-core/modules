process GEOFETCH {
    tag "$geo_accession"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/geofetch:0.12.6--pyh7cba7a3_0':
        'quay.io/biocontainers/geofetch:0.12.6--pyh7cba7a3_0' }"

    input:
    val geo_accession

    output:
    tuple val("${geo_accession}"), path("${geo_accession}/*.CEL.gz"), emit: samples
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    geofetch \\
        -i \\
        $geo_accession \\
        --processed \\
        -g . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        geofetch: \$(geofetch --version|& sed '1!d ; s/geofetch //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    mkdir -p ${geo_accession}
    touch ${geo_accession}/foo.CEL.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        geofetch: \$(geofetch --version|& sed '1!d ; s/geofetch //')
    END_VERSIONS
    """
}
