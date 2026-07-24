process EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eigenstratdatabasetools:1.1.0--hdfd78af_0':
        'quay.io/biocontainers/eigenstratdatabasetools:1.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(geno), path(snp), path(ind)

    output:
    tuple val(meta), path("*.tsv")                                                                                                                      , emit: tsv
    tuple val(meta), path("*.json")                                                                                                                     , emit: json
    tuple val("${task.process}"), val('eigenstratdatabasetools'), eval("eigenstrat_snp_coverage --version 2>&1 | sed 's/^.*eigenstrat_snp_coverage //'"), emit: versions_eigenstratdatabasetools, topic: versions

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
        -j ${prefix}.json \\
        -o ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.json
    """
}
