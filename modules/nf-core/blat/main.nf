process BLAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-blat:472--h9b8f530_0':
        'biocontainers/ucsc-blat:472--h664eb37_1' }"

    input:
    tuple val(meta) , path(query)
    tuple val(meta2), path(subject)

    output:
    tuple val(meta), path("*.psl"), emit: psl
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unzip = query.toString().endsWith(".gz")

    """
    in=$query
    if $unzip
    then
        gunzip -cdf $query > ${prefix}.fasta
        in=${prefix}.fasta
    fi

    blat \\
        $args \\
        $subject \\
        \$in \\
        ${prefix}.psl

    if $unzip
    then
        rm ${prefix}.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blat: \$(echo \$(blat 2>&1) | sed 's/^.*BLAT v. //; s/ fast.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.psl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blat: \$(echo \$(blat 2>&1) | sed 's/^.*BLAT v. //; s/ fast.*\$//')
    END_VERSIONS
    """
}
