process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta), path("*.bed.gz")  , emit: bedgz   , optional: true
    tuple val(meta), path("*.log")     , emit: log     , optional: true
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def reference   = fasta ? "--ref ${fasta}" : ""
    def include_bed = bed ? "--include-bed ${bed}" : ''

    """
    modkit \\
        pileup \\
        $args \\
        --bgzf \\
        --bgzf-threads ${task.cpus} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        $reference \\
        $include_bed \\
        $bam \\
        ${prefix}.tmp

    if test -d ${prefix}.tmp; then
        for file in ${prefix}.tmp/*; do
            if test -f \$file; then
                mv \$file \$(basename \$file)
            fi
        done
    else
        mv ${prefix}.tmp ${prefix}.bed.gz
    fi
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    echo | gzip > ${prefix}.bed.gz
    touch ${prefix}.log
    """
}
