process SAMTOOLS_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input)
    val(interleave)

    output:
    tuple val(meta), path("*_{1,2}.fasta.gz")      , optional:true, emit: fasta
    tuple val(meta), path("*_interleaved.fasta.gz"), optional:true, emit: interleaved
    tuple val(meta), path("*_singleton.fasta.gz")  , optional:true, emit: singleton
    tuple val(meta), path("*_other.fasta.gz")      , optional:true, emit: other
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = ( interleave && ! meta.single_end ) ? "> ${prefix}_interleaved.fasta.gz" :
        meta.single_end ? "-1 ${prefix}_1.fasta.gz -s ${prefix}_singleton.fasta.gz" :
        "-1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz -s ${prefix}_singleton.fasta.gz"
    """
    samtools \\
        fasta \\
        $args \\
        --threads ${task.cpus-1} \\
        -0 ${prefix}_other.fasta.gz \\
        $input \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
