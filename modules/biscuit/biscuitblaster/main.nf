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
    def biscuit_args = task.ext.args ?: ''
    def samblaster_args = task.ext.args2 ?: ''
    def sort_args = task.ext.args3 ?: ''
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    def biscuit_cpus = (int) Math.max(Math.floor(task.cpus/1.05),1)
    def samtools_cpus = task.cpus-biscuit_cpus
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

    biscuit align \\
        -@ $biscuit_cpus \\
        $biscuit_args \\
        $read_group \\
        \$INDEX \\
        $reads | \\
    samblaster \\
        $samblaster_args | \\
    samtools sort \\
        -@ $samtools_cpus \\
        $sort_args \\
        --write-index \\
        -o ${prefix}.bam##idx##${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$(echo \$(biscuit version 2>&1) | sed 's/^.*BISCUIT Version: //; s/Using.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        samblaster: \$(echo \$(samblaster --version 2>&1 | sed 's/^.*samblaster: Version //'))
    END_VERSIONS
    """
}
