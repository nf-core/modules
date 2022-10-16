process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4' : null)
    def container_image = "bowtie2:2.4.4--py39hbb4e92a_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    input:
    path fasta

    output:
    path 'bowtie2'      , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build $args --threads $task.cpus $fasta bowtie2/${fasta.baseName}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
