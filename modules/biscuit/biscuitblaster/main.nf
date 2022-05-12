process BISCUIT_BLASTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113 bioconda::samblaster=0.1.26 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0':
        'quay.io/biocontainers/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0' }"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def biscuit_cpus = (int) Math.max(Math.floor(task.cpus*0.95),1)
    def samtools_cpus = task.cpus-biscuit_cpus
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

    biscuit align \\
        -@ $biscuit_cpus \\
        $args \\
        \$INDEX \\
        $reads | \\
    samblaster \\
        $args2 | \\
    samtools sort \\
        -@ $samtools_cpus \\
        $args3 \\
        --write-index \\
        -o ${prefix}.bam##idx##${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
        samblaster: \$( samblaster --version |& sed 's/^.*samblaster: Version //' )
    END_VERSIONS
    """
}
