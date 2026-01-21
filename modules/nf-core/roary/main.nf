process ROARY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0' :
        'biocontainers/roary:3.13.0--pl526h516909a_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                    , emit: results
    tuple val(meta), path("results/*.aln"), optional: true, emit: aln
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    roary \\
        $args \\
        -p $task.cpus \\
        -f results/ \\
        $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results/
    touch results/core_gene_alignment.aln
    touch results/gene_presence_absence.csv
    touch results/gene_presence_absence.Rtab
    touch results/accessory_binary_genes.fa
    touch results/pan_genome_reference.fa
    touch results/summary_statistics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
    """
}
