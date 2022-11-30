process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::expansionhunter=4.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunter:4.0.2--he785bd8_0' :
        'quay.io/biocontainers/expansionhunter:4.0.2--he785bd8_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path variant_catalog

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gender = (meta.gender == 'male' || meta.gender == 1 || meta.gender == 'XY') ? "male" : "female"
    """
    ExpansionHunter \\
        $args \\
        --reads $bam \\
        --output-prefix $prefix \\
        --reference $fasta \\
        --variant-catalog $variant_catalog \\
        --sex $gender

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$( echo \$(ExpansionHunter --version 2>&1) | sed 's/^.*ExpansionHunter v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$( echo \$(ExpansionHunter --version 2>&1) | sed 's/^.*ExpansionHunter v//')
    END_VERSIONS
    """
}
