process CDHIT_CDHIT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cd-hit=4.8.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit%3A4.8.1--h5b5514e_7':
        'biocontainers/cd-hit:4.8.1--h5b5514e_7' }"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.fasta")    ,emit: fasta
    tuple val(meta), path("*.clstr")    ,emit: clusters
    path "versions.yml"                 ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cd-hit \\
        -i $sequences \\
        -o ${prefix}.fasta \\
        -M $task.memory.mega \\
        -T $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit -h | head -n 1 | sed 's/^.*====== CD-HIT version //;s/ (built on .*) ======//' )
    END_VERSIONS
    """
}
