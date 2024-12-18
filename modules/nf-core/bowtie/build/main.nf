process BOWTIE_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8c0819a9b1f520c49c933e667ae50de2a0730ece4c8b9efe79ac5e403963a9f/data' :
        'community.wave.seqera.io/library/bowtie_samtools:e1a14e1ce4e0170d' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bowtie") , emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bowtie
    bowtie-build --threads $task.cpus $fasta bowtie/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bowtie
    touch bowtie/${prefix}.1.ebwt
    touch bowtie/${prefix}.2.ebwt
    touch bowtie/${prefix}.3.ebwt
    touch bowtie/${prefix}.4.ebwt
    touch bowtie/${prefix}.rev.1.ebwt
    touch bowtie/${prefix}.rev.2.ebwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

}
