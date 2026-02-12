process BOWTIE_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f5ca09fd5aab931d9b87c532c69e0122ce5ff8ec88732f906e12108d48425e9/data' :
        'community.wave.seqera.io/library/bowtie_htslib_samtools:e1e242368ffcb5d3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie') , emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
