process ABACAS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abacas:1.3.1--pl526_0' :
        'biocontainers/abacas:1.3.1--pl526_0' }"

    input:
    tuple val(meta), path(scaffold)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.*"), emit: results
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.abacas"
    """
    abacas.pl \\
        -r ${fasta} \\
        -q ${scaffold} \\
        ${args} \\
        -o ${prefix}

    sort "${prefix}.bin" > "${prefix}.bin.tmp" && mv "${prefix}.bin.tmp" "${prefix}.bin" # Needs to be sorted for consistency
    mv nucmer.delta ${prefix}.nucmer.delta
    mv nucmer.filtered.delta ${prefix}.nucmer.filtered.delta
    mv nucmer.tiling ${prefix}.nucmer.tiling
    mv unused_contigs.out ${prefix}.unused.contigs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abacas: \$(echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}.abacas"
    """
    touch ${prefix}.abacas.MULTIFASTA.fa
    touch ${prefix}.abacas.crunch
    touch ${prefix}.abacas.fasta
    touch ${prefix}.abacas.gaps
    touch ${prefix}.abacas.gaps.tab
    touch ${prefix}.abacas.tab

    touch ${prefix}.nucmer.delta
    touch ${prefix}.nucmer.filtered.delta
    touch ${prefix}.nucmer.tiling
    touch ${prefix}.unused.contigs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abacas: \$(echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//')
    END_VERSIONS
    """
}
