process METHBAT_PROFILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e9081b928a80e8e37f48d53558d39d44ba3b7b05a29055abb3a8e80ca749736/data':
        'community.wave.seqera.io/library/methbat:0.17.0--b493e12136cee7f4' }"

    input:
    tuple val(meta) , path(files)
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("*.tsv"), emit: region_profile
    tuple val(meta), path("*.bed"), emit: asm_bed, optional: true
    tuple val("${task.process}"), val("methbat"), eval("methbat --version 2>&1 | sed 's/.* //'"), emit: versions_methbat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_prefix = files[0].name.toString().replaceAll(/\.((combined)|(hap1)|(hap2))\.bed\.gz(?:\.tbi)?$/, '')
    """
    methbat profile \\
        --input-prefix ${in_prefix} \\
        --input-regions ${regions} \\
        --output-region-profile ${prefix}.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def make_bed = args.contains('--output-asm-bed') ? "touch ${args.split('--output-asm-bed')[1].trim().tokenize(' ')[0]}" : ""
    """
    touch ${prefix}.tsv
    ${make_bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | sed 's/.* //')
    END_VERSIONS
    """
}
