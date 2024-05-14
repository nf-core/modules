process BISCUIT_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::biscuit=1.1.0.20220707 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0':
        'biocontainers/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = (int) Math.max(Math.floor(task.cpus*0.9),1)
    def samtools_cpus = task.cpus-biscuit_cpus
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/\\.bis.amb\$//'`

    biscuit align \\
        $args \\
        -@ $biscuit_cpus \\
        \$INDEX \\
        $reads \\
        | samtools sort $args2 --threads $samtools_cpus --write-index -o ${prefix}.bam##idx##${prefix}.bam.bai -


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
    END_VERSIONS
    """
}
