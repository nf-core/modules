process BCFTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0' :
        'quay.io/biocontainers/bcftools:1.13--h3a49de5_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(fasta)

    output:
    tuple val(meta), path('*.fa'), emit: fasta
    path  "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    cat $fasta | bcftools consensus $vcf $args > ${prefix}.fa
    header=\$(head -n 1 ${prefix}.fa | sed 's/>//g')
    sed -i 's/\${header}/${meta.id}/g' ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
