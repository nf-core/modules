process HISAT3N_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat-3n:0.0.3--b4c98eb79ad7c714' :
        'community.wave.seqera.io/library/hisat-3n:0.0.3--b4b80cb38c483147' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val(meta), path("*.log"), emit: summary
    tuple val("${task.process}"), val('hisat-3n'), eval("hisat-3n --version 2>&1 | head -1 | sed 's/.* //'"), topic: versions, emit: versions_hisat3n

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        INDEX=`find -L ./ -name "*.3n.*.1.ht2" | sed 's/\\.3n\\..*\\.1\\.ht2\$//' | head -1`
        hisat-3n \\
            -x \$INDEX \\
            -U $reads \\
            -S ${prefix}.sam \\
            --summary-file ${prefix}.hisat3n.summary.log \\
            --threads $task.cpus \\
            $args
        """
    } else {
        """
        INDEX=`find -L ./ -name "*.3n.*.1.ht2" | sed 's/\\.3n\\..*\\.1\\.ht2\$//' | head -1`
        hisat-3n \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -S ${prefix}.sam \\
            --summary-file ${prefix}.hisat3n.summary.log \\
            --threads $task.cpus \\
            $args
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sam
    touch ${prefix}.hisat3n.summary.log
    """
}
