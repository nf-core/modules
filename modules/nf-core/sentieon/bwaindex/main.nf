process SENTIEON_BWAINDEX {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80ccb05eb4f1a193a3bd99c4da90f55f74ea6556c25f154e53e1ff5a6caa372d/data' :
        'community.wave.seqera.io/library/sentieon:202503--5e378058d837c58c' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
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
