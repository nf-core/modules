process PHARMCAT_PHENOTYPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.1.1--876b7152770ba008'}"

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
