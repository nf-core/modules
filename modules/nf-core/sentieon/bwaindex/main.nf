process SENTIEON_BWAINDEX {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/sentieon:202308.02--ffce1b7074ce9924' :
        'nf-core/sentieon:202308.02--c641bc397cbf79d5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(bwa), emit: index
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "bwa/${task.ext.prefix}" : "bwa/${fasta.baseName}"
    """
    mkdir bwa

    sentieon \\
        bwa index \\
        $args \\
        -p $prefix \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa

    touch bwa/genome.amb
    touch bwa/genome.ann
    touch bwa/genome.bwt
    touch bwa/genome.pac
    touch bwa/genome.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
