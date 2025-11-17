process METHBAT_PROFILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methbat:0.16.1--h9ee0642_0':
        'biocontainers/methbat:0.16.1--h9ee0642_0' }"

    input:
    tuple val(meta) , path(files)
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("*.tsv"), emit: region_profile
    tuple val(meta), path("*.bed"), emit: asm_bed, optional: true
    path "versions.yml"           , emit: versions

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
