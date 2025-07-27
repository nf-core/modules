process EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eigenstratdatabasetools:1.1.0--hdfd78af_0':
        'biocontainers/eigenstratdatabasetools:1.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(geno), path(snp), path(ind)

    output:
    tuple val(meta), path("*.tsv") , emit: tsv
    tuple val(meta), path("*.json"), emit: json, optional:true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    eigenstrat_snp_coverage \\
        $args \\
        -g ${geno} \\
        -s ${snp} \\
        -i ${ind} \\
        -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eigenstratdatabasetools: \$(echo \$(eigenstrat_snp_coverage --version 2>&1) | sed 's/^.*eigenstrat_snp_coverage //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def json_file = ""
    def matcher = args =~ /-j\s+(\S+)/
    if (matcher.find()) {
        json_file = matcher.group(1)
    }
    def json_cmd = json_file ? "touch ${json_file}" : ""
    """
    touch ${prefix}.tsv
    $json_cmd
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eigenstratdatabasetools: \$(echo \$(eigenstrat_snp_coverage --version 2>&1) | sed 's/^.*eigenstrat_snp_coverage //' ))
    END_VERSIONS
    """
}
