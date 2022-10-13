process BOWTIE_BUILD {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0' : null)
    def container_image = "/bowtie:1.3.0--py38hed8969a_1"
                                               container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
