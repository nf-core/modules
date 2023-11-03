process BOWTIE_BUILD {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "nf-core/modules/bowtie:bowtie-1.3.0_samtools-1.16.1--2c2a8138b83d6d1d"

    input:
    path fasta

    output:
    path 'bowtie'       , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie
    bowtie-build --threads $task.cpus $fasta bowtie/${fasta.baseName}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
