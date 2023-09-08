process MANTA_CONVERTINVERSION {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "bioconda::manta=1.6.0 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:a0332aa38645fbb8969567731ce68cfb7f830ec4-0':
        'biocontainers/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:a0332aa38645fbb8969567731ce68cfb7f830ec4-0' }"

    input:
    tuple val(meta), path(vcf)
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convertInversion.py \$(which samtools) $fasta $vcf | bgzip --threads $task.cpus > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
