process PIRATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.5--hdfd78af_0' :
        'biocontainers/pirate:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${prefix}_results/*")                                   , emit: results
    tuple val(meta), path("${prefix}_results/core_alignment.fasta"), optional: true, emit: aln
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Rename .gff3 to .gff if needed
    find -regex .*\\.gff3\$ | sed 's/3\$//' | xargs -I {} mv {}3 {}

    # Run pirate on all .gff in input directory
    PIRATE \\
        $args \\
        --threads $task.cpus \\
        --input ./ \\
        --output ${prefix}_results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_results
    touch ${prefix}_results/PIRATE.gene_families.ordered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
