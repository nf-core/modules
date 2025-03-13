process BWA_INDEX {
    tag "$fasta"
    // NOTE requires 5.37N memory where N is the size of the database
    // source: https://bio-bwa.sourceforge.net/bwa.shtml#8
    memory { 6.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.18--he4a0461_0' :
        'biocontainers/bwa:0.7.18--he4a0461_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_bwa"), emit: index
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_bwa
    bwa \\
        index \\
        $args \\
        -p ${prefix}_bwa/${fasta.baseName} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir ${prefix}_bwa

    touch ${prefix}_bwa/${fasta.baseName}.amb
    touch ${prefix}_bwa/${fasta.baseName}.ann
    touch ${prefix}_bwa/${fasta.baseName}.bwt
    touch ${prefix}_bwa/${fasta.baseName}.pac
    touch ${prefix}_bwa/${fasta.baseName}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
