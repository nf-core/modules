process BOWTIE_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::bowtie=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hed8969a_1' :
        'biocontainers/bowtie:1.3.0--py38hed8969a_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie') , emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bowtie
    bowtie-build --threads $task.cpus $fasta bowtie/${prefix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
