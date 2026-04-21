process PHARMCAT_PHENOTYPER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7dd711a2b130b55d33e119a346ef8040191bf7834a3c393ed6e29d7d9026d5e/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.2.0--5126bb296d1e59ac'}"

    input:
    tuple val(meta), path(match_json), path(outside_match_tsv)

    output:
    tuple val(meta), path("*.phenotype.json")                                                                   ,   emit: phenotyper_json
    tuple val("${task.process}"), val('pharmcat'), eval("pharmcat --version | cut -f2 -d ' '"), topic: versions ,   emit: versions_pharmcat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}.pharmcat"

    def pheno_input         = match_json        ? " --phenotyper-input ${match_json} "                    : ""
    def outside_pheno_input = outside_match_tsv ? " --phenotyper-outside-call-file ${outside_match_tsv} " : ""

    """
    pharmcat \\
        --base-filename ${prefix} \\
        --output-dir . \\
        -phenotyper \\
        ${pheno_input} \\
        ${outside_pheno_input} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.pharmcat"
    """
    touch ${prefix}.phenotype.json
    """
}
