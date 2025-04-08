process BISCUIT_BLASTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/22f141a92e12f98040054b9a47a51fb4de9c544174bbfa6d251a110f5e26fddc/data':
        'community.wave.seqera.io/library/biscuit_samblaster_samtools:a683d0887e7a91bf' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def biscuit_cpus = [(task.cpus * 0.9) as int, 1].max()
    def samtools_cpus = (task.cpus - biscuit_cpus < 1) ? biscuit_cpus : (task.cpus - biscuit_cpus)
    """
    ln -sf \$(readlink $fasta) $index/$fasta

    biscuit align \\
        -@ $biscuit_cpus \\
        $args \\
        $index/$fasta\\
        $reads \\
        | samblaster $args2 \\
        | samtools sort \\
            -@ $samtools_cpus \\
            $args3 \\
            --write-index \\
            -o ${prefix}.bam -O BAM -

    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
        samblaster: \$( samblaster --version |& sed 's/^.*samblaster: Version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
        samblaster: \$( samblaster --version |& sed 's/^.*samblaster: Version //' )
    END_VERSIONS
    """

}
