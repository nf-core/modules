process BWAMETH_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwameth=0.2.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.2--py_1' :
        'quay.io/biocontainers/bwameth:0.2.2--py_1' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    INDEX=`find -L ${index} -name "*.bwameth.c2t" | sed 's/.bwameth.c2t//'`

    # Modify the timestamps so that bwameth doesn't complain about building the index
    # See https://github.com/nf-core/methylseq/pull/217
    touch -c -- *

    bwameth.py \\
        $args \\
        $read_group \\
        -t $task.cpus \\
        --reference \$INDEX \\
        $reads \\
        | samtools view $args2 -@ $task.cpus -bhS -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(echo \$(bwameth.py --version 2>&1) | cut -f2 -d" ")
    END_VERSIONS
    """
}
