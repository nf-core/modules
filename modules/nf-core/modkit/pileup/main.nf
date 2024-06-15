process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.3.0--h5c23e0d_0':
        'biocontainers/ont-modkit:0.3.0--h5c23e0d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta), path("*.bed")     , emit: bed     , optional: true
    tuple val(meta), path("*.bedgraph"), emit: bedgraph, optional: true
    tuple val(meta), path("*.log")     , emit: log     , optional: true
    path "versions.yml"                , emit: versions

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
        mv ${prefix}.tmp ${prefix}.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/mod_kit //' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.bedgraph
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/mod_kit //' )
    END_VERSIONS
    """
}
