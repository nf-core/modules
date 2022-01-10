process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:cf603b12db30ec91daa04ba45a8ee0f35bbcd1e2-0' :
        'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:cf603b12db30ec91daa04ba45a8ee0f35bbcd1e2-0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa-mem2 \\
        mem \\
        $args \\
        $read_group \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
